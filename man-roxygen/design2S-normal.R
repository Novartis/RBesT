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
#' decision boundary at intermediate points. With the default adaptive
#' scheme (see option \code{RBesT.decision2S_boundary} in
#' \code{\link{RBesT}}) the boundary is evaluated at least once per
#' standard error \eqn{\sigma_2/\sqrt{n_2}} of the mean of sample 2
#' (using the reference scale \eqn{\sigma_2} of prior 2 if
#' \code{family} is given) wherever it is requested and refined further
#' only where it bends, down to the spacing \eqn{\delta_2}. For single
#' component priors with fixed \eqn{\sigma} the boundary is instead an
#' exact line. Operating characteristics and
#' probability of success request the boundary in the regions which
#' hold the probability mass \code{1-eps} of \eqn{y_2}. The decision
#' is represented by a single critical value of \eqn{y_1} for each
#' \eqn{y_2}. If \code{family} is given, then the sampling standard
#' deviation grows without bound for extreme outcomes (e.g.\ response
#' rates close to 0 or 1) such that the posterior reverts to the prior
#' there and the decision can change a second time along \eqn{y_1}.
#' This second change is ignored and results are inaccurate whenever
#' sample 1 has relevant probability mass at such extreme outcomes,
#' which is the case for small expected event counts. The normal
#' approximation is questionable in these settings and the exact
#' designs for the respective endpoint should be used instead. The
#' operating characteristics and the probability of success issue a
#' warning if the affected outcomes have a probability of more than
#' roughly \eqn{0.005}. This probability is estimated on a coarse grid
#' of outcomes and can differ from the actual error by a multiple.
#' @param sigma1 The fixed reference scale of sample 1. If left
#' unspecified, the default reference scale of the prior 1 is assumed.
#' @param sigma2 The fixed reference scale of sample 2. If left
#' unspecified, the default reference scale of the prior 2 is assumed.
#' @param family Optional \code{\link{family}} object specifying a
#' GLM family and link function (e.g.\ \code{binomial()},
#' \code{MASS::negative.binomial(theta)}). When provided, the
#' sampling standard deviation of each sample varies with the
#' respective parameter value via the family's variance function and
#' link. For the Gaussian family \code{sigma1}/\code{sigma2} act as
#' the dispersion parameters and must be supplied; for all other
#' families they must \emph{not} be given (they are determined by the
#' family). Default is \code{NULL} (constant \code{sigma1} and
#' \code{sigma2}).
#' @param offset1 Numeric scalar added to the linear predictor of
#' sample 1 before evaluating the family's variance function.
#' Relevant for Poisson/negative-binomial models with log link where
#' \code{offset1 = log(exposure)}. Default is \code{0}.
#' @param offset2 Numeric scalar added to the linear predictor of
#' sample 2; see \code{offset1}. Defaults to \code{offset1}.
