#' Decision Function for 1 Sample Designs
#'
#' The function sets up a 1 sample decision function with an
#' arbitrary number of conditions.
#'
#' @param pc Vector of critical cumulative probabilities.
#' @param qc Vector of respective critical values. Must match the length of `pc`.
#' @param lower.tail Logical; if `TRUE` (default), probabilities
#' are \eqn{P(X \leq x)}, otherwise, \eqn{P(X > x)}. Either length 1 or same
#' length as `pc`.
#' @param x Two-sided decision function.
#'
#' @details For `lower.tail` being either `TRUE` or `FALSE`,
#' the function creates a one-sided decision function which
#' takes two arguments. The first argument is expected to be a mixture
#' (posterior) distribution. This distribution is tested whether it
#' fulfills all the required threshold conditions specified with the
#' `pc` and `qc` arguments and returns 1 if all conditions
#' are met and 0 otherwise. Hence, for `lower.tail=TRUE`
#' condition \eqn{i} is equivalent to
#'
#' \deqn{P(\theta \leq q_{c,i}) > p_{c,i}}
#'
#' and the decision function is implemented as indicator function on
#' the basis of the heavy-side step function \eqn{H(x)} which is \eqn{0}
#' for \eqn{x \leq 0} and \eqn{1} for \eqn{x > 0}. As all conditions
#' must be met, the final indicator function returns
#'
#' \deqn{\Pi_i H_i(P(\theta \leq q_{c,i}) - p_{c,i} ).}
#'
#' For the case of a boolean vector given to `lower.tail` the
#' direction of each decision aligns respectively, and a two-sided
#' decision function is created.
#'
#' When the second argument is set to `TRUE` a distance metric is
#' returned component-wise per defined condition as
#'
#' \deqn{ D_i = \log(P(\theta < q_{c,i})) - \log(p_{c,i}) .}
#'
#' These indicator functions can be used as input for 1-sample
#' boundary, OC or PoS calculations using [oc1S()] or
#' [pos1S()] .
#'
#' @family design1S
#'
#' @return The function returns a decision function (of class
#' `decision1S_1sided` for one-sided, and of class `decision1S_2sided`
#' for two-sided decisions) which takes two
#' arguments. The first argument is expected to be a mixture
#' (posterior) distribution which is tested if the specified
#' conditions are met. The logical second argument determines if the
#' function acts as an indicator function or if the function returns
#' the distance from the decision boundary for each condition in
#' log-space, i.e. the distance is 0 at the decision boundary,
#' negative for a 0 decision and positive for a 1 decision.
#'
#' For two-sided decision functions, the two components can be
#' extracted with functions [lower()] and [upper()]. The distance
#' as calculated by the decision function is returned as a list with
#' components `lower` and `upper`.
#'
#' @references Neuenschwander B, Rouyrre N, Hollaender H, Zuber E,
#' Branson M. A proof of concept phase II non-inferiority
#' criterion. *Stat. in Med.*. 2011, 30:1618-1627
#'
#' @examples
#'
#' # see Neuenschwander et al., 2011
#'
#' # example is for a time-to-event trial evaluating non-inferiority (NI)
#' # using a normal approximation for the log-hazard ratio
#'
#' # reference scale
#' s <- 2
#' theta_ni <- 0.4
#' theta_a <- 0
#' alpha <- 0.05
#' beta <- 0.2
#' za <- qnorm(1 - alpha)
#' zb <- qnorm(1 - beta)
#' n1 <- round((s * (za + zb) / (theta_ni - theta_a))^2) # n for which design was intended
#' nL <- 233
#' c1 <- theta_ni - za * s / sqrt(n1)
#'
#' # flat prior
#' flat_prior <- mixnorm(c(1, 0, 100), sigma = s)
#'
#' # standard NI design
#' decA <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)
#'
#' # for double criterion with indecision point (mean estimate must be
#' # lower than this)
#' theta_c <- c1
#'
#' # double criterion design
#' # statistical significance (like NI design)
#' dec1 <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)
#' # require mean to be at least as good as theta_c
#' dec2 <- decision1S(0.5, theta_c, lower.tail = TRUE)
#' # combination
#' decComb <- decision1S(c(1 - alpha, 0.5), c(theta_ni, theta_c), lower.tail = TRUE)
#'
#' theta_eval <- c(theta_a, theta_c, theta_ni)
#'
#' # we can display the decision function definition
#' decComb
#'
#' # and use it to decide if a given distribution fulfills all
#' # criterions defined
#' # for the prior
#' decComb(flat_prior)
#' # or for a possible outcome of the trial
#' # here with HR of 0.8 for 40 events
#' decComb(postmix(flat_prior, m = log(0.8), n = 40))
#'
#' # A two-sided decision function can be useful to determine if
#' # certain intermediate (i.e. neither "go" nor "stop") decisions
#' # are to be made based on the posterior distribution.
#' # For example, in the above situation we might have an intermediate
#' # scenario where the trial is significant for non-inferiority but
#' # the mean estimate is in an intermediate range, say between theta_c
#' # theta_f:
#' theta_f <- 0.3
#' decCombIntermediate <- decision1S(
#'   c(1 - alpha, 0.5, 0.8),
#'   c(theta_ni, theta_c, theta_f),
#'   lower.tail = c(TRUE, FALSE, TRUE)
#' )
#' # Not fulfilled for the prior:
#' decCombIntermediate(flat_prior)
#' # But for a hypothetical trial outcome with HR 1.2 and 300 events:
#' decCombIntermediate(postmix(flat_prior, m = log(1.2), n = 300))
#'
#' @export
decision1S <- function(pc = 0.975, qc = 0, lower.tail = TRUE) {
  assert_numeric(pc, lower = 0, upper = 1, any.missing = FALSE, finite = TRUE)
  assert_numeric(qc, len = length(pc), any.missing = FALSE)
  assert_logical(lower.tail, any.missing = FALSE)
  assert_true(length(lower.tail) == 1L || length(lower.tail) == length(pc))
  lower.tail <- scalar_if_same(lower.tail)

  is_two_sided <- length(lower.tail) > 1

  if (is_two_sided) {
    return(create_decision1S_2sided(pc, qc, lower.tail))
  } else {
    return(create_decision1S_1sided(pc, qc, lower.tail))
  }
}

#' @keywords internal
scalar_if_same <- function(x) {
  if (length(x) > 1 && all(x == x[1])) {
    return(x[1])
  }
  x
}

#' Internal Constructor for Atomic 1 Sample One-sided Decision Function
#' @keywords internal
create_decision1S_atomic <- function(pc, qc, lower.tail) {
  lpc <- log(pc)
  fun <- function(mix, dist = FALSE) {
    test <- pmix(mix, qc, lower.tail = lower.tail, log.p = TRUE) - lpc
    if (dist) {
      return(test)
    }
    as.numeric(all(test > 0))
  }
  attr(fun, "pc") <- pc
  attr(fun, "qc") <- qc
  attr(fun, "lower.tail") <- lower.tail

  class(fun) <- c("decision1S_atomic", "function")
  fun
}

#' Internal Constructor for 1 Sample One-sided Decision Function
#' @keywords internal
create_decision1S_1sided <- function(pc, qc, lower.tail) {
  assert_flag(lower.tail)

  atomic_fun <- create_decision1S_atomic(pc, qc, lower.tail)
  attr_name <- if (lower.tail) "lower" else "upper"

  fun <- function(mix, dist = FALSE) {
    test <- atomic_fun(mix, dist)
    if (dist) {
      ret <- stats::setNames(list(test), attr_name)
      return(ret)
    }
    test
  }
  attr(fun, attr_name) <- atomic_fun
  attr(fun, "lower.tail") <- lower.tail

  class(fun) <- c("decision1S", "decision1S_1sided", "function")
  fun
}

#' Internal Constructor for 1 Sample Two-sided Decision Function
#' @keywords internal
create_decision1S_2sided <- function(pc, qc, lower.tail) {
  use_lower <- which(lower.tail)
  use_upper <- which(!lower.tail)
  assert_true(length(use_lower) > 0 && length(use_upper) > 0)

  lower_part <- create_decision1S_atomic(pc[use_lower], qc[use_lower], TRUE)
  upper_part <- create_decision1S_atomic(pc[use_upper], qc[use_upper], FALSE)

  fun <- function(mix, dist = FALSE) {
    dl <- lower_part(mix, dist)
    du <- upper_part(mix, dist)
    if (dist) {
      return(list(lower = dl, upper = du))
    }
    as.numeric(dl && du)
  }
  attr(fun, "lower") <- lower_part
  attr(fun, "upper") <- upper_part

  class(fun) <- c("decision1S", "decision1S_2sided", "function")
  fun
}

#' @rdname decision1S
#' @export
has_lower <- function(x) {
  test_multi_class(x, c("decision1S_2sided", "decision2S_2sided")) ||
    test_multi_class(x, c("decision1S_1sided", "decision2S_1sided")) &&
      attr(x, "lower.tail")
}

#' @rdname decision1S
#' @export
has_upper <- function(x) {
  test_multi_class(x, c("decision1S_2sided", "decision2S_2sided")) ||
    test_multi_class(x, c("decision1S_1sided", "decision2S_1sided")) &&
      !attr(x, "lower.tail")
}

#' @rdname decision1S
#' @export
lower <- function(x) {
  assert_multi_class(x, c("decision1S", "decision2S"))
  attr(x, "lower")
}

#' @rdname decision1S
#' @export
upper <- function(x) {
  assert_multi_class(x, c("decision1S", "decision2S"))
  attr(x, "upper")
}

#' @keywords internal
print_decision1S_atomic <- function(x) {
  qc <- attr(x, "qc")
  pc <- attr(x, "pc")
  low <- attr(x, "lower.tail")
  cmp <- ifelse(low, "<=", ">")
  cat(paste0("P(theta ", cmp, " ", qc, ") > ", pc, "\n"), sep = "")
}

#' @export
print.decision1S_1sided <- function(x, ...) {
  cat("1 sample decision function\n")
  cat("Conditions for acceptance:\n")
  atomic_fun <- if (has_lower(x)) {
    lower(x)
  } else {
    upper(x)
  }
  print_decision1S_atomic(atomic_fun)
  invisible(x)
}

#' @export
print.decision1S_2sided <- function(x, ...) {
  cat("1 sample decision function (two-sided)\n")
  cat("Conditions for acceptance:\n")
  cat("Lower tail conditions:\n")
  print_decision1S_atomic(lower(x))
  cat("Upper tail conditions:\n")
  print_decision1S_atomic(upper(x))
  invisible(x)
}

#' @describeIn decision1S Deprecated old function name. Please use
#' `decision1S` instead.
#' @export
oc1Sdecision <- function(pc = 0.975, qc = 0, lower.tail = TRUE) {
  deprecated("oc1Sdecision", "decision1S")
  return(decision1S(pc, qc, lower.tail))
}
