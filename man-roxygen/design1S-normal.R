#' @describeIn <%= fun %> Applies for the normal model with known
#' standard deviation \eqn{\sigma} and a normal mixture prior for the
#' mean. As a consequence from the assumption of a known standard
#' deviation, the calculation discards sampling uncertainty of the
#' second moment. The function \code{<%= fun %>} has an extra
#' argument \code{eps} (defaults to \eqn{10^{-6}}). The critical value
#' \eqn{y_c} is searched in the region of probability mass
#' \code{1-eps} for \eqn{y}.
#' @param sigma The fixed reference scale. If left unspecified, the
#' default reference scale of the prior is assumed.
