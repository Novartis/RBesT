#' R Bayesian Evidence Synthesis Tools
#'
#' The RBesT tools are designed to support in the derivation of
#' parametric informative priors, asses design characeristics and
#' perform analyses. Supported endpoints include normal, binary and
#' Poisson.
#'
#' For introductory material, please refer to the vignettes which include
#'
#' \itemize{
#' \item Introduction (binary)
#' \item Introduction (normal)
#' \item Customizing RBesT Plots
#' \item Robust MAP, advanced usage
#' }
#'
#' The main function of the package is [gMAP()]. See it's
#' help page for a detailed description of the statistical model.
#'
#'
#' @section Global Options:
#'
#' \tabular{lcl}{
#' Option \tab Default \tab Description \cr
#' `RBesT.MC.warmup` \tab 2000 \tab MCMC warmup iterations \cr
#' `RBesT.MC.iter` \tab 6000 \tab total MCMC iterations \cr
#' `RBesT.MC.chains` \tab 4 \tab MCMC chains\cr
#' `RBesT.MC.thin` \tab 4 \tab MCMC thinning \cr
#' `RBesT.MC.save_warmup` \tab `FALSE` \tab MCMC warmup samples saving \cr
#' `RBesT.MC.control` \tab `list(adapt_delta=0.95,` \tab sets `control` argument for Stan call.\cr
#'  \tab `stepsize=0.01,` \tab `adapt_delta` defaults to `0.99` whenever\cr
#'  \tab `max_treedepth=20)` \tab `RBesT.MC.s2z` is `FALSE`\cr
#' `RBesT.MC.s2z` \tab `TRUE` \tab sum-to-zero parametrization of the group\cr
#' \tab \tab random effects; `FALSE` recovers the legacy\cr
#' \tab \tab sampling scheme\cr
#' `RBesT.MC.ncp` \tab 1 \tab parametrization: 0=CP, 1=NCP, 2=Automatic  \cr
#' `RBesT.MC.init` \tab 1 \tab range of initial uniform \eqn{[-1,1]} is the default  \cr
#' `RBesT.MC.rescale` \tab `TRUE` \tab Automatic rescaling of raw parameters  \cr
#' `RBesT.verbose` \tab `FALSE` \tab requests outputs to be more verbose\cr
#' `RBesT.integrate_args` \tab `list(lower=-Inf,` \tab arguments passed to `integrate` for\cr
#'  \tab `upper=Inf,` \tab adaptive integration of densities (used when\cr
#' \tab `rel.tol=.Machine$double.eps^0.25,` \tab `RBesT.integrate_method` is `"adaptive"`\cr
#' \tab `abs.tol=.Machine$double.eps^0.25,` \tab or for non-`normMix` densities)\cr
#' \tab `subdivisions=1E3)` \tab \cr
#' `RBesT.integrate_prob_eps` \tab `1E-6` \tab probability mass left out from tails if integration needs to be restricted in range \cr
#' `RBesT.integrate_method` \tab `"GQ"` \tab integration method for mixture densities: `"GQ"` (Gaussian quadrature, deterministic) or `"adaptive"` (adaptive Gauss-Kronrod). GQ uses Gauss-Hermite for `normMix`, Gauss-Jacobi for `betaMix`, and Gauss-Laguerre for `gammaMix`. \cr
#' `RBesT.GQ_nodes` \tab `20` \tab starting number of Gaussian quadrature nodes (only used when `RBesT.integrate_method` is `"GQ"`) \cr
#' `RBesT.GQ_rel_tol` \tab `1E-4` \tab relative tolerance target for GQ refinement; the node count is increased until successive estimates agree within `max(RBesT.GQ_abs_tol, RBesT.GQ_rel_tol * |I|)`. Set to a non-finite or non-positive value (e.g. `Inf`) to disable refinement (single evaluation at `RBesT.GQ_nodes`). \cr
#' `RBesT.GQ_abs_tol` \tab `1E-6` \tab absolute tolerance floor for GQ refinement \cr
#' `RBesT.GQ_max_nodes` \tab `240` \tab upper cap on the GQ node count during refinement \cr
#' `RBesT.GQ_node_growth` \tab `2` \tab multiplicative growth factor for the GQ node count between refinement steps \cr
#' `RBesT.GQ_on_nonconvergence` \tab `"adaptive"` \tab behaviour when GQ refinement reaches `RBesT.GQ_max_nodes` without meeting tolerance: `"adaptive"` (fall through to adaptive integration), `"warn"` (warn and return best estimate), `"error"`, or `"silent"` \cr
#' `RBesT.decision2S_boundary` \tab `"adaptive"` \tab tracing scheme for the decision boundary in [decision2S_boundary()] (and hence `oc2S`/`pos2S`): `"adaptive"` (recursive subdivision of the `y2` range with early stop where the boundary is linear for the normal case, and monotone constant-run filling for the binomial/Poisson cases) or `"grid"` (legacy exhaustive/uniform sweep). Boundary values are unchanged (identical for the discrete families, within solver tolerance for the normal case); `"adaptive"` needs far fewer root solves. \cr
#' }
#'
#' @section Version History:
#'
#' See `NEWS.md` file.
#'
#' @references
#' \insertRef{rstan}{RBesT}
#'
#' @useDynLib RBesT, .registration = TRUE
#'
# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @import abind
#' @import assertthat
#' @import checkmate
#' @import Formula
#' @import ggplot2
#' @import methods
#' @import mvtnorm
#' @import Rcpp
#' @import rstantools
#' @import stats
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom lifecycle deprecated
#' @importFrom matrixStats rowLogSumExps colLogSumExps colSums2 rowMins rowRanks logSumExp
#' @importFrom RcppParallel RcppParallelLibs CxxFlags
#' @importFrom rlang .data
#' @importFrom rstan sampling extract get_sampler_params summary
#' @importFrom posterior as_draws_array as_draws_rvars as_draws_matrix as_draws_list as_draws_df subset_draws
#' @importFrom Rdpack reprompt
#' @importFrom utils capture.output modifyList
#' @export nsamples
## usethis namespace: end
"_PACKAGE"
