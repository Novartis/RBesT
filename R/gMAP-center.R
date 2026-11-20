#' Construct working Gaussian data for partial centering
#'
#' Approximates each supported likelihood on its linear-predictor scale and
#' aggregates the resulting information by random-effects group.
#'
#' @param family Character string naming the likelihood family.
#' @param y,y.se Gaussian observations and their standard errors.
#' @param r,n Binomial event counts and sample sizes.
#' @param count,log.offset Poisson counts and log offsets.
#' @param group.index Integer group index for each observation.
#' @param n.groups Number of random-effects groups.
#'
#' @return A list containing the working response, standard error, observation
#'   information, and information summed by group.
#' @noRd
.gmap_working_data <- function(family, y, y.se, r, n, count, log.offset,
                               group.index, n.groups) {
  if (family == "gaussian") {
    response <- y
    se <- y.se
  } else if (family == "binomial") {
    p <- (r + 0.5) / (n + 1)
    response <- qlogis(p)
    se <- 1 / sqrt(n * p * (1 - p))
  } else if (family == "poisson") {
    response <- log(count + 0.5) - log.offset
    se <- 1 / sqrt(count + 0.5)
  } else {
    stop("Unsupported likelihood family for partial centering: ", family)
  }

  response <- as.numeric(response)
  se <- as.numeric(se)
  information <- 1 / se^2
  group.information <- numeric(n.groups)
  grouped.information <- rowsum(information, group.index, reorder = FALSE)
  group.information[as.integer(rownames(grouped.information))] <-
    grouped.information
  list(response = response, se = se, information = information,
       group.information = group.information)
}

#' Evaluate a heterogeneity-prior quantile function
#'
#' @param p Numeric vector of probabilities.
#' @param dist Character string naming the heterogeneity prior.
#' @param prior Numeric vector containing the prior parameters.
#'
#' @return The prior quantiles at `p`, restricted to non-negative
#'   heterogeneity where required.
#' @noRd
.gmap_tau_quantile <- function(p, dist, prior) {
  p1 <- prior[1]
  p2 <- prior[2]
  switch(
    dist,
    HalfNormal = qnorm((1 + p) / 2, sd = p2),
    TruncNormal = {
      lower.cdf <- pnorm(0, mean = p1, sd = p2)
      qnorm(lower.cdf + p * (1 - lower.cdf), mean = p1, sd = p2)
    },
    Uniform = {
      lower <- max(0, p1)
      lower + p * (p2 - lower)
    },
    Gamma = qgamma(p, shape = p1, rate = p2),
    InvGamma = 1 / qgamma(1 - p, shape = p1, rate = p2),
    LogNormal = qlnorm(p, meanlog = p1, sdlog = p2),
    TruncCauchy = {
      lower.cdf <- pcauchy(0, location = p1, scale = p2)
      qcauchy(lower.cdf + p * (1 - lower.cdf),
              location = p1, scale = p2)
    },
    Exp = qexp(p, rate = p1),
    stop("Unsupported tau prior for partial centering: ", dist)
  )
}

#' Construct a tensor-product quadrature rule for heterogeneity
#'
#' @param n Number of Gauss-Legendre nodes in each dimension.
#' @param quantile.functions List of heterogeneity-prior quantile functions,
#'   one per active stratum.
#'
#' @return A list containing heterogeneity values, their corresponding prior
#'   probabilities, and tensor-product quadrature weights.
#' @noRd
.gmap_tensor_rule <- function(n, quantile.functions) {
  rule <- .get_GQ_rule(n, kind = "legendre")
  # Map the Legendre rule from [-1, 1] to prior probabilities in [0, 1].
  probability <- (rule$nodes + 1) / 2
  one.weight <- rule$weights / 2
  dimension <- length(quantile.functions)
  # Every row selects one node per active heterogeneity stratum.
  index <- do.call(
    expand.grid,
    c(rep(list(seq_len(n)), dimension), KEEP.OUT.ATTRS = FALSE)
  )
  index <- as.matrix(index)

  tau <- vapply(
    seq_len(dimension),
    function(j) quantile.functions[[j]](probability[index[, j]]),
    numeric(nrow(index))
  )
  if (dimension == 1L) {
    tau <- matrix(tau, ncol = 1L)
  }
  node.probability <- vapply(
    seq_len(dimension),
    function(j) probability[index[, j]],
    numeric(nrow(index))
  )
  if (dimension == 1L) {
    node.probability <- matrix(node.probability, ncol = 1L)
  }
  weight <- apply(index, 1L, function(i) prod(one.weight[i]))
  list(tau = tau, probability = node.probability, weight = weight)
}

#' Build the working-model marginal log likelihood
#'
#' Integrates the Gaussian fixed effects analytically and returns a function
#' that evaluates the resulting marginal log likelihood at heterogeneity
#' quadrature nodes.
#'
#' @param work Working Gaussian data from `.gmap_working_data()`.
#' @param X Fixed-effects design matrix.
#' @param beta.prior Matrix of fixed-effect prior means and standard deviations.
#' @param group.index Integer group index for each observation.
#' @param group.stratum Integer heterogeneity stratum for each group.
#' @param active.strata Strata represented in the quadrature rule.
#'
#' @return A function accepting a matrix of heterogeneity nodes and returning
#'   their marginal log likelihoods.
#' @noRd
.gmap_working_marginal <- function(work, X, beta.prior, group.index,
                                   group.stratum, active.strata) {
  residual <- as.numeric(work$response - X %*% beta.prior[, 1])
  scaled.X <- sweep(X, 2L, beta.prior[, 2], "*")
  V <- cbind(residual, scaled.X)
  inv.d <- 1 / work$se^2
  Vinv <- sweep(V, 1L, inv.d, "*")
  C0 <- crossprod(V, Vinv)

  n.groups <- length(group.stratum)
  A <- matrix(0, nrow = n.groups, ncol = ncol(V))
  grouped.A <- rowsum(Vinv, group.index, reorder = FALSE)
  A[as.integer(rownames(grouped.A)), ] <- grouped.A
  group.precision <- numeric(n.groups)
  grouped.precision <- rowsum(inv.d, group.index, reorder = FALSE)
  group.precision[as.integer(rownames(grouped.precision))] <-
    grouped.precision
  stratum.column <- match(group.stratum, active.strata)
  active.group <- !is.na(stratum.column)
  log.det.base <- sum(log(work$se^2))
  n.beta <- ncol(X)

  function(tau) {
    apply(tau, 1L, function(tau.node) {
      tau.group <- numeric(n.groups)
      tau.group[active.group] <- tau.node[stratum.column[active.group]]
      # Evaluate the marginal MVN using the Woodbury matrix identity and the
      # matrix determinant lemma. With D = diag(work$se^2) and Z the
      # observation-by-group incidence matrix, the direct equivalent is:
      # mvtnorm::dmvnorm(
      #   residual,
      #   sigma = D + Z %*% diag(tau.group^2) %*% t(Z) +
      #     tcrossprod(scaled.X),
      #   log = TRUE
      # )
      # It is not used because it factorizes a full observation-sized
      # covariance matrix at every quadrature node; the result below differs
      # only by the tau-independent MVN normalizing constant.
      shrinkage <- tau.group^2 / (1 + tau.group^2 * group.precision)
      # Apply the grouped random-effect covariance update to the crossproduct.
      C <- C0 - crossprod(A, sweep(A, 1L, shrinkage, "*"))
      M <- C[-1L, -1L, drop = FALSE] + diag(n.beta)
      chol.M <- tryCatch(chol(M), error = function(e) NULL)
      if (is.null(chol.M)) {
        return(-Inf)
      }
      b <- C[-1L, 1L]
      solved <- backsolve(chol.M, b, transpose = TRUE)
      quadratic <- C[1L, 1L] - sum(solved^2)
      log.det <- log.det.base +
        sum(log1p(tau.group^2 * group.precision)) +
        2 * sum(log(diag(chol.M)))
      -0.5 * (log.det + max(0, quadratic))
    })
  }
}

#' Extract a marginal posterior median from a tensor rule
#'
#' Uses interpolation in prior-probability space for the main estimate and a
#' weighted discrete quantile as an independent numerical guard.
#'
#' @param x Heterogeneity values for one stratum.
#' @param probability Prior probabilities corresponding to `x`.
#' @param mass Normalized joint posterior mass at every tensor node.
#' @param rule.weight Tensor-product quadrature weights.
#'
#' @return A named numeric vector containing the interpolated `main` median and
#'   the discrete `guard` median.
#' @noRd
.gmap_marginal_median <- function(x, probability, mass, rule.weight) {
  keep <- is.finite(x) & is.finite(probability) &
    is.finite(mass) & mass > 0 &
    is.finite(rule.weight) & rule.weight > 0
  x <- x[keep]
  probability <- probability[keep]
  mass <- mass[keep]
  rule.weight <- rule.weight[keep]
  if (length(x) < 2L) {
    stop("insufficient finite quadrature support")
  }

  order.x <- order(x)
  x <- x[order.x]
  probability <- probability[order.x]
  mass <- mass[order.x]
  rule.weight <- rule.weight[order.x]
  # Collapse repeated marginal nodes contributed by the other dimensions.
  unique.index <- match(x, unique(x))
  x <- unique(x)
  weighted.probability <- as.numeric(rowsum(
    probability * rule.weight, unique.index, reorder = FALSE
  ))
  mass <- as.numeric(rowsum(mass, unique.index, reorder = FALSE))
  rule.weight <- as.numeric(rowsum(rule.weight, unique.index,
                                   reorder = FALSE))
  probability <- weighted.probability / rule.weight
  # Removing the integration weight recovers density on the probability scale.
  log.density <- log(mass) - log(rule.weight)
  # Exclude negligible tails before interpolating across widely spaced quantiles.
  retain <- log.density >= max(log.density) - 50
  x.main <- x[retain]
  probability.main <- probability[retain]
  log.density <- log.density[retain]
  if (length(x.main) < 2L || diff(range(x.main)) <= 0) {
    stop("insufficient posterior support after density trimming")
  }

  grid <- seq(
    min(probability.main), max(probability.main), length.out = 4001L
  )
  density <- exp(approx(probability.main, log.density, xout = grid,
                        rule = 2, ties = "ordered")$y -
                   max(log.density))
  increment <- diff(grid) * (density[-length(density)] + density[-1L]) / 2
  cdf <- c(0, cumsum(increment))
  if (!is.finite(cdf[length(cdf)]) || cdf[length(cdf)] <= 0) {
    stop("invalid interpolated posterior normalizer")
  }
  cdf <- cdf / cdf[length(cdf)]
  median.probability <- approx(
    cdf, grid, xout = 0.5, ties = "ordered"
  )$y
  main <- approx(
    probability.main, x.main, xout = median.probability,
    rule = 2, ties = "ordered"
  )$y

  # This discrete estimate guards against interpolation artifacts.
  normalized.mass <- mass / sum(mass)
  mass.midpoint <- cumsum(normalized.mass) - normalized.mass / 2
  guard <- approx(
    mass.midpoint, x, xout = 0.5, rule = 2, ties = "ordered"
  )$y
  c(main = main, guard = guard)
}

#' Approximate posterior heterogeneity medians by quadrature
#'
#' @param n Number of Gauss-Legendre nodes in each dimension.
#' @param quantile.functions List of heterogeneity-prior quantile functions.
#' @param log.marginal Function evaluating the marginal log likelihood.
#'
#' @return A list containing main and guard medians for each stratum and the
#'   number of likelihood evaluations.
#' @noRd
.gmap_tau_posterior <- function(n, quantile.functions, log.marginal) {
  rule <- .gmap_tensor_rule(n, quantile.functions)
  log.likelihood <- log.marginal(rule$tau)
  if (!any(is.finite(log.likelihood))) {
    stop("all quadrature likelihood evaluations are non-finite")
  }
  log.mass <- log(rule$weight) + log.likelihood
  log.normalizer <- matrixStats::logSumExp(log.mass)
  if (!is.finite(log.normalizer)) {
    stop("invalid quadrature normalizer")
  }
  mass <- exp(log.mass - log.normalizer)
  median <- vapply(
    seq_len(ncol(rule$tau)),
    function(j) .gmap_marginal_median(
      rule$tau[, j], rule$probability[, j], mass, rule$weight
    ),
    numeric(2L)
  )
  if (ncol(rule$tau) == 1L) {
    median <- matrix(median, nrow = 2L,
                     dimnames = list(c("main", "guard"), NULL))
  }
  list(main = median["main", ], guard = median["guard", ],
       evaluations = nrow(rule$tau), rule = rule, mass = mass)
}

#' Build the conditional-Gaussian evaluator at a fixed heterogeneity node
#'
#' Extends the Woodbury reduction already used by `.gmap_working_marginal()`
#' to also return, at every heterogeneity node, the conditional posterior
#' mean/covariance of `beta` (marginal over the group effects) and the
#' conditional posterior mean/variance of every group effect (the pure,
#' zero-mean random deviation Stan calls `re`; marginal over both the
#' residual sampling variance and the `beta` uncertainty, via the law of
#' total variance). This reduction is part of the production quadrature
#' so no second, independent implementation of the Woodbury step exists.
#'
#' @inheritParams .gmap_working_marginal
#' @return A function of one heterogeneity node (a numeric vector, one entry
#'   per active stratum) returning `NULL` on a non-positive-definite fit and
#'   otherwise a list with `beta.mean`, `beta.cov`, `group.mean`,
#'   `group.var`.
#' @noRd
.gmap_conditional_posterior <- function(work, X, beta.prior, group.index,
                                        group.stratum, active.strata,
                                        t.df = Inf) {
  residual <- as.numeric(work$response - X %*% beta.prior[, 1])
  scaled.X <- sweep(X, 2L, beta.prior[, 2], "*")
  V <- cbind(residual, scaled.X)
  inv.d <- 1 / work$se^2
  Vinv <- sweep(V, 1L, inv.d, "*")
  C0 <- crossprod(V, Vinv)

  n.groups <- length(group.stratum)
  n.beta <- ncol(X)
  A <- matrix(0, nrow = n.groups, ncol = ncol(V))
  grouped.A <- rowsum(Vinv, group.index, reorder = FALSE)
  A[as.integer(rownames(grouped.A)), ] <- grouped.A
  group.precision <- numeric(n.groups)
  grouped.precision <- rowsum(inv.d, group.index, reorder = FALSE)
  group.precision[as.integer(rownames(grouped.precision))] <-
    grouped.precision
  stratum.column <- match(group.stratum, active.strata)
  active.group <- !is.na(stratum.column)
  # Student-t random effects: the same local-curvature surrogate already used
  # for the centering fraction (no exact state-independent Fisher fraction
  # exists for Student-t).
  curvature <- if (is.finite(t.df)) t.df / (t.df + 1) else 1

  function(tau.node) {
    tau.group <- numeric(n.groups)
    tau.group[active.group] <- tau.node[stratum.column[active.group]]
    tau2 <- tau.group^2 * curvature
    shrinkage <- tau2 / (1 + tau2 * group.precision)
    shrinkage[group.precision == 0] <- tau2[group.precision == 0]
    C <- C0 - crossprod(A, sweep(A, 1L, shrinkage, "*"))
    M <- C[-1L, -1L, drop = FALSE] + diag(n.beta)
    chol.M <- tryCatch(chol(M), error = function(e) NULL)
    if (is.null(chol.M)) {
      return(NULL)
    }
    b <- C[-1L, 1L]
    z <- backsolve(chol.M, b, transpose = TRUE)
    mu.std <- backsolve(chol.M, z)
    cov.std <- chol2inv(chol.M)

    beta.mean <- beta.prior[, 1] + beta.prior[, 2] * mu.std
    beta.cov <- outer(beta.prior[, 2], beta.prior[, 2]) * cov.std

    # group_j | beta, tau, y ~ N(fraction_j * r_j, shrinkage_j), with
    # r_j = ybar_j - xbar_j'beta linear in beta_std; marginalizing beta adds
    # fraction_j^2 * Var(r_j) to the conditional variance (total variance).
    fraction <- group.precision * shrinkage
    coef <- A[, -1, drop = FALSE]
    has.info <- group.precision > 0
    coef[has.info, ] <- coef[has.info, , drop = FALSE] / group.precision[has.info]
    Er <- numeric(n.groups)
    Er[has.info] <- (A[has.info, 1] - as.numeric(
      A[has.info, -1, drop = FALSE] %*% mu.std
    )) / group.precision[has.info]
    Vr <- rowSums((coef %*% cov.std) * coef)
    group.mean <- fraction * Er
    group.var <- shrinkage + fraction^2 * Vr

    list(beta.mean = beta.mean, beta.cov = beta.cov,
         group.mean = group.mean, group.var = group.var)
  }
}

#' Aggregate conditional posteriors over a quadrature rule
#'
#' Applies the law of total expectation/variance to combine
#' `.gmap_conditional_posterior()` evaluations across a set of heterogeneity
#' nodes and normalized posterior masses, for both a full tensor rule and the
#' degenerate one-node case (`tau.dist == "Fixed"`, mass = 1).
#'
#' @param tau.nodes Matrix of heterogeneity nodes, one row per node.
#' @param mass Normalized posterior mass at each node (sums to one).
#' @param n.groups,mX Dimensions of the group and fixed-effect vectors.
#'
#' @return A list with `ok = FALSE` if no node produced a positive-definite
#'   fit, otherwise `ok = TRUE` and `beta.location`, `beta.scale`,
#'   `group.location`, `group.scale`.
#' @noRd
.gmap_conditional_aggregate <- function(conditional, tau.nodes, mass,
                                        n.groups, mX) {
  nodes <- lapply(seq_len(nrow(tau.nodes)), function(i) {
    conditional(tau.nodes[i, ])
  })
  valid <- !vapply(nodes, is.null, logical(1))
  if (!any(valid)) {
    return(list(ok = FALSE))
  }
  nodes <- nodes[valid]
  w <- mass[valid] / sum(mass[valid])
  n.node <- length(nodes)

  beta.means <- matrix(vapply(nodes, function(x) x$beta.mean, numeric(mX)),
                       nrow = n.node, ncol = mX, byrow = TRUE)
  beta.vars <- matrix(vapply(nodes, function(x) diag(x$beta.cov), numeric(mX)),
                      nrow = n.node, ncol = mX, byrow = TRUE)
  beta.location <- as.numeric(w %*% beta.means)
  beta.var <- as.numeric(w %*% beta.vars) +
    as.numeric(w %*% (sweep(beta.means, 2, beta.location)^2))

  group.means <- matrix(
    vapply(nodes, function(x) x$group.mean, numeric(n.groups)),
    nrow = n.node, ncol = n.groups, byrow = TRUE
  )
  group.vars <- matrix(
    vapply(nodes, function(x) x$group.var, numeric(n.groups)),
    nrow = n.node, ncol = n.groups, byrow = TRUE
  )
  group.location <- as.numeric(w %*% group.means)
  group.var <- as.numeric(w %*% group.vars) +
    as.numeric(w %*% (sweep(group.means, 2, group.location)^2))

  list(ok = TRUE, beta.location = beta.location,
       beta.scale = sqrt(pmax(beta.var, 0)),
       group.location = group.location,
       group.scale = sqrt(pmax(group.var, 0)))
}

#' Choose partial-centering fractions for a gMAP model
#'
#' Estimates posterior heterogeneity under a working Gaussian model, checks
#' quadrature stability, and maps heterogeneity to group-specific centering
#' fractions.
#'
#' @param family Character string naming the likelihood family.
#' @param y,y.se Gaussian observations and their standard errors.
#' @param r,n Binomial event counts and sample sizes.
#' @param count,log.offset Poisson counts and log offsets.
#' @param X Fixed-effects design matrix.
#' @param beta.prior Matrix of fixed-effect prior means and standard deviations.
#' @param group.index Integer group index for each observation.
#' @param group.stratum Integer heterogeneity stratum for each group.
#' @param tau.dist Character string naming the heterogeneity prior.
#' @param tau.prior Matrix of heterogeneity-prior parameters by stratum.
#' @param prior_PD Logical indicating a prior-predictive model.
#' @param t.df Degrees of freedom for Student-t random effects.
#'
#' @return A list reporting whether the approximation succeeded (`ok`), any
#'   failure `reason`, total quadrature `evaluations`, the partial-centering
#'   fractions (`center`) and heterogeneity medians (`tau`) already consumed
#'   by `RBesT.MC.ncp` 2/3, and the additional sampler-scaling summaries
#'   consumed by `R/gMAP.R`: per-stratum `log.tau.location`/`log.tau.scale`,
#'   `beta.location`/`beta.scale`, and `group.location`/`group.scale`. Fine
#'   and coarse tensor rules are evaluated once each and reused for every one
#'   of these outputs -- no output triggers a second, independent quadrature
#'   pass.
#' @noRd
.gmap_quadrature_approximation <- function(family, y, y.se, r, n, count,
                                           log.offset, X, beta.prior,
                                           group.index, group.stratum,
                                           tau.dist, tau.prior, prior_PD,
                                           t.df = Inf) {
  n.groups <- length(group.stratum)
  n.strata <- nrow(tau.prior)
  mX <- ncol(X)
  empty.result <- function(reason = NULL, group.information = numeric(n.groups)) {
    list(ok = is.null(reason), reason = reason, evaluations = 0L,
         center = numeric(n.groups), tau = numeric(n.strata),
         log.tau.location = rep(NA_real_, n.strata),
         log.tau.scale = rep(NA_real_, n.strata),
         beta.location = rep(NA_real_, mX), beta.scale = rep(NA_real_, mX),
         group.location = numeric(n.groups),
         group.scale = rep(NA_real_, n.groups),
         group.information = group.information)
  }
  if (prior_PD) {
    return(empty.result())
  }

  work <- .gmap_working_data(
    family, y, y.se, r, n, count, log.offset, group.index, n.groups
  )
  observed.strata <- sort(unique(group.stratum[work$group.information > 0]))
  if (any(vapply(observed.strata, function(s) {
    sum(work$group.information > 0 & group.stratum == s) < 2L
  }, logical(1)))) {
    return(empty.result("a tau stratum has fewer than two informative groups", work$group.information))
  }

  conditional <- .gmap_conditional_posterior(
    work, X, beta.prior, group.index, group.stratum, observed.strata, t.df
  )

  if (tau.dist == "Fixed") {
    tau.median <- tau.prior[, 1]
    posterior.guard <- tau.median
    evaluations <- 0L
    log.tau.location <- rep(NA_real_, n.strata)
    log.tau.scale <- rep(NA_real_, n.strata)
    fixed.tau <- matrix(tau.median[observed.strata], nrow = 1L)
    aggregate <- .gmap_conditional_aggregate(
      conditional, fixed.tau, 1, n.groups, mX
    )
    if (!aggregate$ok) {
      return(empty.result("non-positive-definite conditional posterior at fixed tau", work$group.information))
    }
  } else {
    if (length(observed.strata) > 3L) {
      return(empty.result("more than three observed non-fixed tau strata", work$group.information))
    }
    quantile.functions <- lapply(observed.strata, function(s) {
      force(s)
      function(p) .gmap_tau_quantile(p, tau.dist, tau.prior[s, ])
    })
    log.marginal <- .gmap_working_marginal(
      work, X, beta.prior, group.index, group.stratum, observed.strata
    )
    # Keep the tensor rule tractable as the number of strata increases.
    schedule <- switch(as.character(length(observed.strata)),
                       "1" = c(fine = 48L, coarse = 32L),
                       "2" = c(fine = 32L, coarse = 24L),
                       "3" = c(fine = 24L, coarse = 20L))
    fine <- tryCatch(
      .gmap_tau_posterior(schedule["fine"], quantile.functions, log.marginal),
      error = function(e) e
    )
    coarse <- tryCatch(
      .gmap_tau_posterior(schedule["coarse"], quantile.functions, log.marginal),
      error = function(e) e
    )
    if (inherits(fine, "error") || inherits(coarse, "error")) {
      reason <- if (inherits(fine, "error")) conditionMessage(fine) else
        conditionMessage(coarse)
      return(empty.result(reason, work$group.information))
    }
    if (any(!is.finite(c(fine$main, fine$guard, coarse$main)))) {
      return(empty.result("non-finite posterior tau median", work$group.information))
    }
    # Agreement across rule sizes detects quadrature that has not converged.
    if (any(abs(log(fine$main / coarse$main)) > 0.1)) {
      return(empty.result("fine and coarse posterior tau medians disagree", work$group.information))
    }
    tau.median <- numeric(n.strata)
    posterior.guard <- numeric(n.strata)
    tau.median[observed.strata] <- fine$main
    posterior.guard[observed.strata] <- fine$guard
    evaluations <- fine$evaluations + coarse$evaluations

    # log(tau) location/scale: reuse the fine rule's nodes and mass already
    # computed for the main/guard medians above -- no extra evaluation of
    # log.marginal().
    log.tau <- log(fine$rule$tau)
    log.tau.mean <- as.numeric(fine$mass %*% log.tau)
    log.tau.var <- as.numeric(fine$mass %*% (sweep(log.tau, 2, log.tau.mean)^2))
    log.tau.location <- rep(NA_real_, n.strata)
    log.tau.scale <- rep(NA_real_, n.strata)
    log.tau.location[observed.strata] <- log.tau.mean
    log.tau.scale[observed.strata] <- sqrt(pmax(log.tau.var, 0))

    # beta/group location-scale: aggregate over the same fine rule.
    aggregate <- .gmap_conditional_aggregate(
      conditional, fine$rule$tau, fine$mass, n.groups, mX
    )
    if (!aggregate$ok) {
      return(empty.result("no positive-definite quadrature node", work$group.information))
    }
  }

  curvature <- if (is.finite(t.df)) t.df / (t.df + 1) else 1
  make.center <- function(tau) {
    variance <- tau[group.stratum]^2 * curvature
    information <- work$group.information
    fraction <- variance * information / (1 + variance * information)
    fraction[information == 0 | tau[group.stratum] == 0] <- 0
    fraction
  }
  center <- make.center(tau.median)
  guard.center <- make.center(posterior.guard)
  if (any(!is.finite(center)) || any(center < 0 | center > 1)) {
    return(empty.result("invalid partial-centering fractions", work$group.information))
  }
  if (any(abs(center - guard.center) > 0.01)) {
    return(empty.result("posterior median extractors disagree", work$group.information))
  }

  if (tau.dist != "Fixed") {
    coarse.tau <- numeric(n.strata)
    coarse.tau[observed.strata] <- coarse$main
    if (any(abs(center - make.center(coarse.tau)) > 0.01)) {
      return(empty.result("fine and coarse centering fractions disagree", work$group.information))
    }
  }
  list(ok = TRUE, center = center, tau = tau.median, reason = NULL,
       evaluations = evaluations,
       log.tau.location = log.tau.location, log.tau.scale = log.tau.scale,
       beta.location = aggregate$beta.location,
       beta.scale = aggregate$beta.scale,
       group.location = aggregate$group.location,
       group.scale = aggregate$group.scale,
       group.information = work$group.information)
}
