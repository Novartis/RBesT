#' @param prior1 Prior for sample 1.
#' @param prior2 Prior for sample 2.
#' @param n1,n2 Sample size of the respective samples. Sample size \code{n1} must be greater than 0 while sample size \code{n2} must be greater or equal to 0. 
#' @param decision Two-sample decision function to use; see \code{\link{decision2S}}.
#' @param ... Optional arguments.
#' @param eps Support of random variables are determined as the
#' interval covering \code{1-eps} probability mass. Defaults to
#' \eqn{10^{-6}}.
#' @param Ngrid Determines density of discretization grid on which
#' decision function is evaluated (see below for more details).
