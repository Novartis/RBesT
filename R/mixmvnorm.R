#' @name mixmvnorm
#'
#' @title Multivariate Normal Mixture Density
#'
#' @description The multivariate normal mixture density and auxiliary
#'     functions.
#'
#' @param ... List of mixture components.
#' @param sigma Reference covariance.
#' @param param Determines how the parameters in the list are
#'     interpreted. See details.
#' @param object,x Multivariate normal mixture object.
#'
#' @details Each entry in the `...` argument list is a numeric
#'     vector defining one component of the mixture multivariate
#'     normal distribution. The first entry of the component defining
#'     vector is the weight of the mixture component followed by the
#'     vector of means in each dimension and finally a specification
#'     of the covariance matrix, which depends on the chosen
#'     parametrization. The covariance matrix is expected to be given
#'     as numeric vector in a column-major format, which is standard
#'     conversion applied to matrices by the vector concatenation
#'     function [base::c()]. Please refer to the examples
#'     section below.
#'
#' Each component defining vector can be specified in different ways
#' as determined by the `param` option:
#'
#' \describe{
#' \item{ms}{Mean vector and covariance matrix `s`. Default.}
#' \item{mn}{Mean vector and number of observations. `n` determines
#' the covariance for each component via the relation \eqn{\Sigma/n}
#' with \eqn{\Sigma} being the known reference covariance.}
#' \item{msr}{Mean vector, standard deviations and correlations in
#' column-major format (corresponds to order when printing multi-variate
#' normal mixtures).}
#' }
#'
#' The reference covariance \eqn{\Sigma} is the known covariance in
#' the normal-normal model (observation covariance). The function
#' `sigma` can be used to query the reference covariance and may
#' also be used to assign a new reference covariance, see examples
#' below. In case `sigma` is not specified, the user has to
#' supply `sigma` as argument to functions which require a
#' reference covariance.
#'
#' @family mixdist
#'
#' @return Returns a multivariate normal mixture with the specified
#'     mixture components.
#'
#' @examples
#'
#' # default mean & covariance parametrization
#' S <- diag(c(1, 2)) %*% matrix(c(1, 0.5, 0.5, 1), 2, 2) %*% diag(c(1, 2))
#' mvnm1 <- mixmvnorm(
#'   rob = c(0.2, c(0, 0), diag(c(2, 2)^2)),
#'   inf = c(0.8, c(0.5, 1), S / 4), sigma = S
#' )
#'
#' print(mvnm1)
#' summary(mvnm1)
#'
#' set.seed(657846)
#' mixSamp1 <- rmix(mvnm1, 500)
#' colMeans(mixSamp1)
#'
#' # alternative mean, sd and correlation parametrization
#' mvnm1_alt <- mixmvnorm(
#'   rob = c(0.2, c(0, 0), c(2, 2), 0.0),
#'   inf = c(0.8, c(0.5, 1), c(1, 2) / 2, 0.5),
#'   sigma = msr2mvnorm(s = c(1, 2), r = 0.5, unlist = FALSE)$s,
#'   param = "msr"
#' )
#'
#' print(mvnm1_alt)
#'
NULL

#' @rdname mixmvnorm
#' @export
mixmvnorm <- function(..., sigma, param = c("ms", "mn", "msr")) {
  param <- match.arg(param)
  ## length of first mean vector determines dimension
  mix <- mixdist3(...)
  dim_labels <- rownames(mix)
  Nc <- ncol(mix)
  n <- colnames(mix)
  if (param == "ms") {
    ## mean vector & covariance parametrization
    l <- nrow(mix)
    p <- (sqrt(1 + 4 * (l - 1)) - 1) / 2
    assert_integerish(p, lower = 1, any.missing = FALSE, len = 1)
    ## in this case we expect c(weight, mean, as.numeric(cov))
    mix <- do.call(
      mixdist3,
      lapply(
        1:Nc,
        function(co)
          c(
            mix[1, co],
            mvnorm(
              mix[2:(p + 1), co],
              matrix(mix[(1 + p + 1):(1 + p + p^2), co], p, p)
            )
          )
      )
    )
  }
  if (param == "mn") {
    ## mean vector & number of observations
    l <- nrow(mix)
    p <- l - 2
    assert_integerish(p, lower = 1, any.missing = FALSE, len = 1)
    assert_matrix(sigma, any.missing = FALSE, nrows = p, ncols = p)
    mix <- do.call(
      mixdist3,
      lapply(
        1:Nc,
        function(co) {
          assert_numeric(
            mix[l, co],
            lower = 0,
            finite = TRUE,
            any.missing = FALSE
          )
          c(mix[1, co], mvnorm(mix[2:(p + 1), co], sigma / mix[l, co]))
        }
      )
    )
  }
  if (param == "msr") {
    ## mean vector, sds & correlations in column-major ordering
    l <- nrow(mix)
    p <- -3 / 2 + sqrt(9 / 4 - 2 * (1 - l))
    assert_integerish(p, lower = 1, any.missing = FALSE, len = 1)
    ## in this case we expect c(weight, mean, sds, rho)
    mix <- do.call(
      mixdist3,
      lapply(
        1:Nc,
        function(co) {
          mvnc <- mix[, co]
          names(mvnc)[-1] <- mvnorm_label(mix[-1, co])
          mvnc
        }
      )
    )
  }
  assert_numeric(
    mix[1, ],
    lower = 0,
    upper = 1,
    finite = TRUE,
    any.missing = FALSE,
    .var.name = "weights"
  )
  colnames(mix) <- n
  p <- mvnormdim(mix[-1, 1])
  if (is.null(dim_labels)) {
    dim_labels <- as.character(1:p)
  } else {
    dim_labels <- dim_labels[2:(p + 1)]
  }
  rownames(mix) <- c("w", mvnorm_label(mix[-1, 1], dim_labels))
  if (!missing(sigma)) {
    assert_matrix(sigma, any.missing = FALSE, nrows = p, ncols = p)
    colnames(sigma) <- rownames(sigma) <- dim_labels
    attr(mix, "sigma") <- sigma
  }
  class(mix) <- c("mvnormMix", "mix")
  likelihood(mix) <- "mvnormal"
  mix
}

#' @rdname mixmvnorm
#' @param m Mean vector.
#' @param s Standard deviation vector.
#' @param r Vector of correlations in column-major format of the lower
#'   triangle of the correlation matrix.
#' @param unlist Logical. Controls whether the result is a flattened
#'   vector (`TRUE`) or a list with mean `m` and covariance `s`
#'   (`FALSE`). Defaults to `TRUE`.
#' @export
msr2mvnorm <- function(m, s, r, unlist = TRUE) {
  if (unlist) {
    checkmate::assert_numeric(
      m,
      finite = TRUE,
      any.missing = FALSE,
      min.len = 1
    )
    d <- length(m)
  } else {
    if (missing(m)) {
      d <- length(s)
    } else {
      d <- length(m)
    }
  }
  checkmate::assert_numeric(
    s,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    len = d
  )
  n_r <- d * (d - 1) / 2
  Rho <- diag(1, d, d)
  if (n_r == 0) {
    assert_that(
      missing(r) || is.null(r) || length(r) == 0,
      msg = "Assertion on 'r' failed: r must be NULL or missing if dimension is one."
    )
  }
  if (n_r > 0) {
    checkmate::assert_numeric(
      r,
      lower = -1,
      upper = 1,
      finite = TRUE,
      any.missing = FALSE,
      len = n_r
    )
    Rho[lower.tri(Rho)] <- r
  }
  Rho[upper.tri(Rho)] <- Rho[lower.tri(Rho)]
  cov <- diag(s, d, d) %*% Rho %*% diag(s, d, d)
  if (unlist) {
    return(c(m, cov))
  }
  if (!missing(m)) {
    return(list(m = m, s = cov))
  }
  list(s = cov)
}


#' @keywords internal
mvnorm <- function(mean, sigma) {
  ## TODO: far more checks!!! Allow to pass in directly a cholesky factor
  assert_numeric(mean, finite = TRUE, any.missing = FALSE)
  p <- length(mean)
  assert_matrix(sigma, any.missing = FALSE, nrows = p, ncols = p)

  # Compute standard deviations
  s <- sqrt(diag(sigma))

  # Handle zero standard deviations
  zero_sd <- s == 0
  if (any(zero_sd)) {
    rho <- diag(1, nrow = p, ncol = p) # Initialize correlation matrix with 0
    if (!all(zero_sd)) {
      rho[!zero_sd, !zero_sd] <- cov2cor(sigma[
        !zero_sd,
        !zero_sd,
        drop = FALSE
      ]) # Compute correlations for non-zero SDs
    }
  } else {
    rho <- cov2cor(sigma) # Standard case
  }

  # Construct the result
  mvn <- c(mean, s, rho[lower.tri(rho)])
  names(mvn) <- mvnorm_label(mvn)
  mvn
}


#' @keywords internal
mvnormdim <- function(mvn) {
  p <- (-3 + sqrt(9 + 8 * length(mvn))) / 2
  assert_integerish(p, lower = 1, any.missing = FALSE, len = 1)
  as.integer(p)
}

#' @keywords internal
mvnorm_dim_labels <- function(mvn) {
  p <- mvnormdim(mvn)
  n <- names(mvn)
  if (is.null(n)) {
    return(as.character(1:p))
  }
  return(gsub("\\]$", "", gsub("^[^\\[]*\\[", "", n[1:p])))
}

#' @keywords internal
mvnorm_label <- function(mvn, dim_labels) {
  p <- mvnormdim(mvn)
  if (missing(dim_labels)) {
    dim_labels <- mvnorm_dim_labels(mvn)
  }
  if (p > 1) {
    Rho_labs_idx <- outer(
      dim_labels,
      dim_labels,
      paste,
      sep = ","
    )[lower.tri(diag(p))]
    lab <- c(
      paste0("m[", dim_labels, "]"),
      paste0("s[", dim_labels, "]"),
      paste0("rho[", Rho_labs_idx, "]")
    )
  } else {
    lab <- c(paste0("m[", dim_labels, "]"), paste0("s[", dim_labels, "]"))
  }
  lab
}

#' @keywords internal
mvnormsigma <- function(mvn) {
  p <- mvnormdim(mvn)
  n <- length(mvn)
  s <- mvn[(p + 1):(2 * p)]
  Rho <- diag(p)
  Rho[lower.tri(Rho)] <- mvn[(2 * p + 1):n]
  Rho[upper.tri(Rho)] <- t(Rho)[upper.tri(Rho)]
  diag(s, nrow = p) %*% Rho %*% diag(s, nrow = p)
}


#' @rdname mixmvnorm
#' @method print mvnormMix
#' @export
print.mvnormMix <- function(x, ...) {
  cat("Multivariate normal mixture\n")
  cat(paste0("Outcome dimension: ", mvnormdim(x[-1, 1]), "\n"))
  if (!is.null(sigma(x))) {
    cat("Reference covariance:\n")
    print(sigma(x), ...)
  }
  NextMethod()
}

#' @rdname mixmvnorm
#' @method summary mvnormMix
#' @export
summary.mvnormMix <- function(object, ...) {
  w <- object[1, ]
  p <- mvnormdim(object[-1, 1])
  m <- object[2:(p + 1), , drop = FALSE]
  Nc <- ncol(object)
  mmix <- rowSums(sweep(m, 2, w, "*"))
  ## Cov(x,x) = E(x x') - E(x) E(x')
  ## E(x x') = Sigma + m m' (see matrix cookbook eq 377)
  if (Nc == 1) {
    S <- w[1] * mvnormsigma(object[-1, 1])
  } else {
    S <- -1 * tcrossprod(mmix)
    for (i in 1:Nc) {
      S <- S +
        w[i] *
          (mvnormsigma(object[-1, i]) +
            tcrossprod(unname(m[, i, drop = FALSE])))
    }
  }
  rownames(S) <- colnames(S) <- names(mmix) <- mvnorm_dim_labels(object[-1, 1])
  list(mean = mmix, cov = S)
}

#' @rdname mixmvnorm
#' @method sigma mvnormMix
#' @export
#' @export sigma
sigma.mvnormMix <- function(object, ...) {
  attr(object, "sigma")
}

#' @keywords internal
is_mixmv <- function(mix) {
  inherits(mix, "mvnormMix")
}

#' @keywords internal
mvnorm_extract_dim <- function(mix, sub) {
  Nc <- ncol(mix)
  sub_comp <- list()
  p <- mvnormdim(mix[-1, 1])
  assert_numeric(sub, lower = 1, upper = p, any.missing = FALSE)
  for (i in seq_len(Nc)) {
    sub_comp[[i]] <- c(
      mix["w", i],
      mix[1 + sub, i],
      mvnormsigma(mix[-1, i])[sub, sub, drop = FALSE]
    )
  }
  if (!is.null(sigma(mix))) {
    sub_comp$sigma <- sigma(mix)
  }
  do.call(mixmvnorm, sub_comp)
}
