#' Decision Function for 2 Sample Designs
#'
#' The function sets up a 2 sample one-sided decision function with an
#' arbitrary number of conditions on the difference distribution.
#'
#' @param pc Vector of critical cumulative probabilities of the
#' difference distribution.
#' @param qc Vector of respective critical values of the difference
#' distribution. Must match the length of `pc`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities
#' are \eqn{P(X \leq x)}, otherwise, \eqn{P(X > x)}.
#' @param link Enables application of a link function prior to
#' evaluating the difference distribution. Can take one of the values
#' `identity` (default), `logit` or `log`.
#'
#' @details This function creates a decision function (of class `decision2S`
#' for one-sided, and of class `decision2S_2sided` for two-sided decisions)
#' on the basis of the difference distribution in a 2 sample situation. To
#' support double criterion designs, see *Neuenschwander et al.,
#' 2010*, an arbitrary number of criterions can be given. The decision
#' function demands that the probability mass below the critical value
#' `qc` of the difference \eqn{\theta_1 - \theta_2} is at least
#' `pc`. Hence, for `lower.tail=TRUE` condition \eqn{i} is
#' equivalent to
#'
#' \deqn{P(\theta_1 - \theta_2 \leq q_{c,i}) > p_{c,i}}
#'
#' and the decision function is implemented as indicator function
#' using the heavy-side step function \eqn{H(x)} which is \eqn{0} for
#' \eqn{x \leq 0} and \eqn{1} for \eqn{x > 0}. As all conditions must
#' be met, the final indicator function returns
#'
#' \deqn{\Pi_i H_i(P(\theta_1 - \theta_2 \leq q_{c,i}) - p_{c,i} ),}
#'
#' which is \eqn{1} if all conditions are met and \eqn{0}
#' otherwise. For `lower.tail=FALSE` differences must be greater
#' than the given quantiles `qc`.
#'
#' For the case of a boolen vector given to `lower.tail` the
#' direction of each decision aligns respectively, and a two-sided
#' decision function is created.
#'
#' Note that whenever a `link` other than `identity` is
#' requested, then the underlying densities are first transformed
#' using the link function and then the probabilties for the
#' differences are calculated in the transformed space. Hence, for a
#' binary endpoint the default `identity` link will calculate
#' risk differences, the `logit` link will lead to decisions
#' based on the differences in `logit`s corresponding to a
#' criterion based on the log-odds. The `log` link will evaluate
#' ratios instead of absolute differences which could be useful for a
#' binary endpoint or counting rates. The respective critical
#' quantiles `qc` must be given on the transformed scale.
#'
#' @return The function returns a decision function which takes three
#' arguments. The first and second argument are expected to be mixture
#' (posterior) distributions from which the difference distribution is
#' formed and all conditions are tested. The third argument determines
#' if the function acts as an indicator function or if the function
#' returns the distance from the decision boundary for each condition
#' in log-space. That is, the distance is 0 at the decision boundary,
#' negative for a 0 decision and positive for a 1 decision.
#'
#' For two-sided decision functions, the two components can be
#' extracted with functions [lower()] and [upper()]. The distance
#' as calculated by the decision function is returned as a list with
#' components `lower` and `upper`.
#'
#' @references Gsponer T, Gerber F, Bornkamp B, Ohlssen D,
#' Vandemeulebroecke M, Schmidli H.A practical guide to Bayesian group
#' sequential designs. *Pharm. Stat.*. 2014; 13: 71-80
#'
#' @family design2S
#'
#' @examples
#'
#' # see Gsponer et al., 2010
#' priorT <- mixnorm(c(1, 0, 0.001), sigma = 88, param = "mn")
#' priorP <- mixnorm(c(1, -49, 20), sigma = 88, param = "mn")
#' # the success criteria is for delta which are larger than some
#' # threshold value which is why we set lower.tail=FALSE
#' successCrit <- decision2S(c(0.95, 0.5), c(0, 50), FALSE)
#' # the futility criterion acts in the opposite direction
#' futilityCrit <- decision2S(c(0.90), c(40), TRUE)
#'
#' print(successCrit)
#' print(futilityCrit)
#'
#' # consider decision for specific outcomes
#' postP_interim <- postmix(priorP, n = 10, m = -50)
#' postT_interim <- postmix(priorT, n = 20, m = -80)
#' futilityCrit(postP_interim, postT_interim)
#' successCrit(postP_interim, postT_interim)
#'
#' # Binary endpoint with double criterion decision on log-odds scale
#' # 95% certain positive difference and an odds ratio of 2 at least
#' decL2 <- decision2S(c(0.95, 0.5), c(0, log(2)), lower.tail = FALSE, link = "logit")
#' # 95% certain positive difference and an odds ratio of 3 at least
#' decL3 <- decision2S(c(0.95, 0.5), c(0, log(3)), lower.tail = FALSE, link = "logit")
#'
#' # data scenario
#' post1 <- postmix(mixbeta(c(1, 1, 1)), n = 40, r = 10)
#' post2 <- postmix(mixbeta(c(1, 1, 1)), n = 40, r = 18)
#'
#' # positive outcome and a median odds ratio of at least 2 ...
#' decL2(post2, post1)
#' # ... but not more than 3
#' decL3(post2, post1)
#'
#' @export
decision2S <- function(
  pc = 0.975,
  qc = 0,
  lower.tail = TRUE,
  link = c("identity", "logit", "log")
) {
  assert_numeric(pc)
  assert_numeric(qc, len = length(pc))
  assert_logical(lower.tail)
  assert_true(length(lower.tail) == 1L || length(lower.tail) == length(pc))
  lower.tail <- scalar_if_same(lower.tail)
  link <- match.arg(link)

  is_two_sided <- length(lower.tail) > 1

  if (is_two_sided) {
    create_decision2S_2sided(pc, qc, lower.tail, link)
  } else {
    create_decision2S_1sided(pc, qc, lower.tail, link)
  }
}

#' Internal Constructor for 2 Sample One-sided Decision Function
#' @keywords internal
create_decision2S_1sided <- function(pc, qc, lower.tail, link) {
  lpc <- log(pc)
  dlink_obj <- link_map[[link]]
  fun <- function(mix1, mix2, dist = FALSE) {
    dlink(mix1) <- dlink_obj
    dlink(mix2) <- dlink_obj
    ## Note that for normal mixture densities we can expedite the
    ## calculation of the convolution dramatically, i.e. the
    ## convolution is done analytically exact
    test <- if (inherits(mix1, "normMix")) {
      pmix(mixnormdiff(mix1, mix2), qc, lower.tail = lower.tail, log.p = TRUE) -
        lpc
    } else {
      log(pmax(
        pmixdiff(mix1, mix2, qc, lower.tail = lower.tail),
        .Machine$double.eps
      )) -
        lpc
    }
    if (dist) {
      return(test)
    }
    as.numeric(all(test > 0))
  }
  attr(fun, "pc") <- pc
  attr(fun, "qc") <- qc
  attr(fun, "link") <- link
  attr(fun, "lower.tail") <- scalar_if_same(lower.tail)

  class(fun) <- c("decision2S", "function")
  fun
}

#' Internal Constructor for 2 Sample Two-sided Decision Function
#' @keywords internal
create_decision2S_2sided <- function(pc, qc, lower.tail, link) {
  use_lower <- which(lower.tail)
  use_upper <- which(!lower.tail)
  assert_true(length(use_lower) > 0 && length(use_upper) > 0)

  lower_part <- create_decision2S_1sided(
    pc[use_lower],
    qc[use_lower],
    TRUE,
    link
  )
  upper_part <- create_decision2S_1sided(
    pc[use_upper],
    qc[use_upper],
    FALSE,
    link
  )

  fun <- function(mix1, mix2, dist = FALSE) {
    dl <- lower_part(mix1, mix2, dist)
    du <- upper_part(mix1, mix2, dist)
    if (dist) {
      return(list(lower = dl, upper = du))
    }
    as.numeric(all(dl > 0) && all(du > 0))
  }
  attr(fun, "lower") <- lower_part
  attr(fun, "upper") <- upper_part
  attr(fun, "link") <- link

  class(fun) <- c("decision2S_2sided", "function")
  fun
}

#' @keywords internal
print_decision2S_1sided <- function(x) {
  qc <- attr(x, "qc")
  pc <- attr(x, "pc")
  low <- attr(x, "lower.tail")
  cmp <- ifelse(low, "<=", ">")
  cat(paste0("P(theta1 - theta2 ", cmp, " ", qc, ") > ", pc, "\n"), sep = "")
}

#' @export
print.decision2S <- function(x, ...) {
  cat("2 sample decision function\n")

  cat("Conditions for acceptance:\n")
  print_decision2S_1sided(x)

  link <- attr(x, "link")
  cat("Link:", link, "\n")

  invisible(x)
}

#' @export
print.decision2S_2sided <- function(x, ...) {
  cat("2 sample two-sided decision function\n")

  cat("Lower side conditions for acceptance:\n")
  print_decision2S_1sided(lower(x))
  cat("Upper side conditions for acceptance:\n")
  print_decision2S_1sided(upper(x))

  link <- attr(x, "link")
  cat("Link:", link, "\n")
  invisible(x)
}

#' @describeIn decision2S Deprecated old function name. Please use
#' `decision2S` instead.
#' @export
oc2Sdecision <- function(
  pc = 0.975,
  qc = 0,
  lower.tail = TRUE,
  link = c("identity", "logit", "log")
) {
  deprecated("oc2Sdecision", "decision2S")
  return(decision2S(pc, qc, lower.tail, link))
}
