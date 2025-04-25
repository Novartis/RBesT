#' Operating Characteristics for 2 Sample Design
#'
#' The `oc2S` function defines a 2 sample design (priors, sample
#' sizes & decision function) for the calculation of operating
#' characeristics. A function is returned which calculates the
#' calculates the frequency at which the decision function is
#' evaluated to 1 when assuming known parameters.
#'
#' @template args-boundary2S
#'
#' @details The `oc2S` function defines a 2 sample design and
#' returns a function which calculates its operating
#' characteristics. This is the frequency with which the decision
#' function is evaluated to 1 under the assumption of a given true
#' distribution of the data defined by the known parameter
#' \eqn{\theta_1} and \eqn{\theta_2}. The 2 sample design is defined
#' by the priors, the sample sizes and the decision function,
#' \eqn{D(y_1,y_2)}. These uniquely define the decision boundary , see
#' [decision2S_boundary()].
#'
#' Calling the `oc2S` function calculates the decision boundary
#' \eqn{D_1(y_2)} (see [decision2S_boundary()]) and returns
#' a function which can be used to calculate the desired frequency
#' which is evaluated as
#'
#' \deqn{ \int f_2(y_2|\theta_2) F_1(D_1(y_2)|\theta_1) dy_2. }
#'
#' See below for examples and specifics for the supported mixture
#' priors.
#'
#' @return Returns a function which when called with two arguments
#' `theta1` and `theta2` will return the frequencies at
#' which the decision function is evaluated to 1 whenever the data is
#' distributed according to the known parameter values in each
#' sample. Note that the returned function takes vector arguments.
#'
#' @references Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B.
#' Robust meta-analytic-predictive priors in clinical trials with historical control information.
#' *Biometrics* 2014;70(4):1023-1032.
#'
#' @family design2S
#'
#' @examples
#'
#' # example from Schmidli et al., 2014
#' dec <- decision2S(0.975, 0, lower.tail = FALSE)
#'
#' prior_inf <- mixbeta(c(1, 4, 16))
#' prior_rob <- robustify(prior_inf, weight = 0.2, mean = 0.5)
#' prior_uni <- mixbeta(c(1, 1, 1))
#'
#' N <- 40
#' N_ctl <- N - 20
#'
#' # compare designs with different priors
#' design_uni <- oc2S(prior_uni, prior_uni, N, N_ctl, dec)
#' design_inf <- oc2S(prior_uni, prior_inf, N, N_ctl, dec)
#' design_rob <- oc2S(prior_uni, prior_rob, N, N_ctl, dec)
#'
#' # type I error
#' curve(design_inf(x, x), 0, 1)
#' curve(design_uni(x, x), lty = 2, add = TRUE)
#' curve(design_rob(x, x), lty = 3, add = TRUE)
#'
#' # power
#' curve(design_inf(0.2 + x, 0.2), 0, 0.5)
#' curve(design_uni(0.2 + x, 0.2), lty = 2, add = TRUE)
#' curve(design_rob(0.2 + x, 0.2), lty = 3, add = TRUE)
#'
#' @export
oc2S <- function(prior1, prior2, n1, n2, decision, ...) UseMethod("oc2S")
#' @export
oc2S.default <- function(prior1, prior2, n1, n2, decision, ...)
  "Unknown density"

#' @templateVar fun oc2S
#' @template design2S-binomial
#' @export
oc2S.betaMix <- function(prior1, prior2, n1, n2, decision, eps, ...) {
  if (missing(eps) & ((n1 + 1) * (n2 + 1) > 1e7)) {
    warning("Large sample space. Consider setting eps=1e-6.")
  }

  crit_y1 <- decision2S_boundary(prior1, prior2, n1, n2, decision, eps)

  lower.tail <- attr(decision, "lower.tail")

  approx_method <- !missing(eps)

  design_fun <- function(theta1, theta2, y2) {
    ## other-wise we calculate the frequencies at which the
    ## decision is 1 (probability mass with decision==1)

    ## in case n2==0, then theta2 is irrelevant
    if (n2 == 0 & missing(theta2)) {
      theta2 <- 0.5
    }

    if (!missing(y2)) {
      deprecated("Use of y2 argument", "decision2S_boundary")
      return(crit_y1(y2, lim1 = c(0, n1)))
    }

    assert_numeric(theta1, lower = 0, upper = 1, finite = TRUE)
    assert_numeric(theta2, lower = 0, upper = 1, finite = TRUE)

    T <- try(data.frame(theta1 = theta1, theta2 = theta2, row.names = NULL))
    if (inherits(T, "try-error")) {
      stop("theta1 and theta2 need to be of same size")
    }

    ## now get the decision boundary in the needed range
    lim1 <- c(0, n1)
    lim2 <- c(0, n2)

    ## in case we use the approximate method, we restrict the
    ## evaluation of the decision function range
    if (approx_method) {
      lim1[1] <- qbinom(eps / 2, n1, min(theta1))
      lim1[2] <- qbinom(1 - eps / 2, n1, max(theta1))
      lim2[1] <- qbinom(eps / 2, n2, min(theta2))
      lim2[2] <- qbinom(1 - eps / 2, n2, max(theta2))
    }

    ## for each 0:n1 of the possible outcomes, calculate the
    ## probability mass past the boundary (log space) weighted with
    ## the density as the value for 1 occures (due to theta1)
    boundary <- crit_y1(lim2[1]:lim2[2], lim1 = c(0, n1))
    res <- matrix(-Inf, nrow = diff(lim2) + 1, ncol = nrow(T))

    for (i in lim2[1]:lim2[2]) {
      y2ind <- i - lim2[1] + 1
      if (boundary[y2ind] == -1) {
        ## decision was always 0
        res[y2ind, ] <- -Inf
      } else if (boundary[y2ind] == n1 + 1) {
        ## decision was always 1
        res[y2ind, ] <- 0
      } else {
        ## calculate for all requested theta1 the probability mass
        ## past (or before) the boundary
        res[y2ind, ] <- pbinom(
          boundary[y2ind],
          n1,
          T$theta1,
          lower.tail = lower.tail,
          log.p = TRUE
        )
      }
      ## finally weight with the density according to the occurence
      ## of i due to theta2; the pmax avoids -Inf in a case of Prob==0
      ## res[y2ind,] <- res[y2ind,] + pmax(dbinom(i, n2, T$theta2, log=TRUE), -700)
      res[y2ind, ] <- res[y2ind, ] + dbinom(i, n2, T$theta2, log = TRUE)
    }
    ## exp(log_colSum_exp(res))
    exp(matrixStats::colLogSumExps(res))
  }
  design_fun
}

#' @templateVar fun oc2S
#' @template design2S-normal
#' @export
oc2S.normMix <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps = 1e-6,
  Ngrid = 10,
  ...
) {
  ## distributions of the means of the data generating distributions
  ## for now we assume that the underlying standard deviation
  ## matches the respective reference scales
  if (missing(sigma1)) {
    sigma1 <- RBesT::sigma(prior1)
    message("Using default prior 1 reference scale ", sigma1)
  }
  if (missing(sigma2)) {
    sigma2 <- RBesT::sigma(prior2)
    message("Using default prior 2 reference scale ", sigma2)
  }
  assert_number(sigma1, lower = 0)
  assert_number(sigma2, lower = 0)

  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2

  crit_y1 <- decision2S_boundary(
    prior1,
    prior2,
    n1,
    n2,
    decision,
    sigma1,
    sigma2,
    eps,
    Ngrid
  )

  lower.tail <- attr(decision, "lower.tail")

  sem1 <- sigma1 / sqrt(n1)
  sem2 <- sigma2 / sqrt(n2)

  if (n2 == 0) sem2 <- sigma(prior2) / sqrt(1E-1)

  ## change the reference scale of the prior such that the prior
  ## represents the distribution of the respective means
  mean_prior1 <- prior1
  sigma(mean_prior1) <- sem1
  ## mean_prior2 <- prior2
  ## sigma(mean_prior2) <- sem2

  freq <- function(theta1, theta2) {
    lim1 <- qnorm(c(eps / 2, 1 - eps / 2), theta1, sem1)
    lim2 <- qnorm(c(eps / 2, 1 - eps / 2), theta2, sem2)
    if (n2 == 0) {
      return(pnorm(crit_y1(theta2), theta1, sem1, lower.tail = lower.tail))
    } else {
      return(integrate_density_log(
        function(x)
          pnorm(
            crit_y1(x, lim1 = lim1),
            theta1,
            sem1,
            lower.tail = lower.tail,
            log.p = TRUE
          ),
        mixnorm(c(1, theta2, sem2), sigma = sem2),
        logit(eps / 2),
        logit(1 - eps / 2)
      ))
    }
  }

  Vfreq <- Vectorize(freq)

  design_fun <- function(theta1, theta2, y2) {
    if (!missing(y2)) {
      theta2 <- y2
    }

    ## in case n2==0, then theta2 is irrelevant
    if (n2 == 0 & missing(theta2)) {
      theta2 <- theta1
    }

    lim2 <- c(
      qnorm(p = eps / 2, mean = min(theta2), sd = sem2),
      qnorm(p = 1 - eps / 2, mean = max(theta2), sd = sem2)
    )

    if (!missing(y2)) {
      deprecated("Use of y2 argument", "decision2S_boundary")
      return(crit_y1(y2))
    }

    ## ensure that boundary is calculated for the full range
    ## needed
    lim1 <- c(
      qnorm(eps / 2, min(theta1), sem1),
      qnorm(1 - eps / 2, max(theta1), sem1)
    )
    lim2 <- c(
      qnorm(eps / 2, min(theta2), sem2),
      qnorm(1 - eps / 2, max(theta2), sem2)
    )

    ## call boundary function to cache all results for all
    ## requested computations
    crit_y1(lim2, lim1 = lim1)

    T <- try(data.frame(theta1 = theta1, theta2 = theta2, row.names = NULL))
    if (inherits(T, "try-error")) {
      stop("theta1 and theta2 need to be of same size")
    }

    do.call(Vfreq, T)
  }

  design_fun
}

#' @templateVar fun oc2S
#' @template design2S-poisson
#' @export
oc2S.gammaMix <- function(prior1, prior2, n1, n2, decision, eps = 1e-6, ...) {
  assert_that(likelihood(prior1) == "poisson")
  assert_that(likelihood(prior2) == "poisson")

  crit_y1 <- decision2S_boundary(prior1, prior2, n1, n2, decision, eps)

  lower.tail <- attr(decision, "lower.tail")

  freq <- function(theta1, theta2) {
    lambda1 <- theta1 * n1
    lambda2 <- theta2 * n2
    lim1 <- qpois(c(eps / 2, 1 - eps / 2), lambda1)
    grid <- seq(qpois(eps / 2, lambda2), qpois(1 - eps / 2, lambda2))

    exp(matrixStats::logSumExp(
      dpois(grid, lambda2, log = TRUE) +
        ppois(
          crit_y1(grid, lim1 = lim1),
          lambda1,
          lower.tail = lower.tail,
          log.p = TRUE
        )
    ))
  }

  Vfreq <- Vectorize(freq)

  design_fun <- function(theta1, theta2, y2) {
    ## in case n2==0, then theta2 is irrelevant
    if (n2 == 0 & missing(theta2)) {
      theta2 <- theta1
    }

    if (!missing(y2)) {
      if (missing(theta1)) {
        theta1 <- summary(prior1)["mean"]
      }
      lambda2 <- y2
    } else {
      lambda2 <- theta2 * n2
    }

    lambda1 <- theta1 * n1

    lim1 <- c(0, 0)
    lim1[1] <- qpois(eps / 2, min(lambda1))
    lim1[2] <- qpois(1 - eps / 2, max(lambda1))

    if (n2 == 0) {
      lim2 <- c(0, 0)
    } else {
      lim2 <- c(0, 0)
      lim2[1] <- qpois(eps / 2, min(lambda2))
      lim2[2] <- qpois(1 - eps / 2, max(lambda2))
    }

    ## force lower limit of lim1 to be 0 such that we will get and
    ## answer in most cases; performance wise it should be ok as
    ## we run a O(log(N)) search
    lim1[1] <- 0

    ## ensure that the boundaries are cached
    crit_y1(lim2, lim1 = lim1)

    if (!missing(y2)) {
      deprecated("Use of y2 argument", "decision2S_boundary")
      return(crit_y1(y2, lim1 = lim1))
    }

    T <- try(data.frame(theta1 = theta1, theta2 = theta2, row.names = NULL))
    if (inherits(T, "try-error")) {
      stop("theta1 and theta2 need to be of same size")
    }
    do.call(Vfreq, T)
  }
  design_fun
}
