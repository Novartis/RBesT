#' @describeIn <%= fun %> Applies for the normal model with known
#' standard deviation \eqn{\sigma} and normal mixture priors for the
#' means. As a consequence from the assumption of a known standard
#' deviation, the calculation discards sampling uncertainty of the
#' second moment. The function has two extra arguments (with
#' defaults): \code{eps} (\eqn{10^{-6}}) and \code{Ngrid} (10). The
#' decision boundary is searched in the region of probability mass
#' \code{1-eps}, respectively for \eqn{y_1} and \eqn{y_2}. The
#' continuous decision function is evaluated at a discrete grid, which
#' is determined by a spacing with \eqn{\delta_2 =
#' \sigma_2/\sqrt{N_{grid}}}. Once the decision boundary is evaluated
#' at the discrete steps, a spline is used to inter-polate the
#' decision boundary at intermediate points.
#' @param sigma1 The fixed reference scale of sample 1. If left
#' unspecified, the default reference scale of the prior 1 is assumed.
#' @param sigma2 The fixed reference scale of sample 2. If left
#' unspecified, the default reference scale of the prior 2 is assumed.
