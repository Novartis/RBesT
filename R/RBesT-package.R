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
#' `RBesT.MC.control` \tab `list(adapt_delta=0.99,` \tab sets `control` argument for Stan call\cr
#'  \tab `stepsize=0.01,` \tab \cr
#'  \tab `max_treedepth=20)` \tab \cr
#' `RBesT.MC.ncp` \tab 1 \tab parametrization: 0=CP, 1=NCP, 2=Automatic  \cr
#' `RBesT.MC.init` \tab 1 \tab range of initial uniform \eqn{[-1,1]} is the default  \cr
#' `RBesT.MC.rescale` \tab `TRUE` \tab Automatic rescaling of raw parameters  \cr
#' `RBesT.verbose` \tab `FALSE` \tab requests outputs to be more verbose\cr
#' `RBesT.integrate_args` \tab `list(lower=-Inf,` \tab arguments passed to `integrate` for\cr
#'  \tab `upper=Inf,` \tab intergation of densities\cr
#' \tab `rel.tol=.Machine$double.eps^0.25,` \tab \cr
#' \tab `abs.tol=.Machine$double.eps^0.25,` \tab \cr
#' \tab `subdivisions=1E3)` \tab \cr
#' `RBesT.integrate_prob_eps` \tab `1E-6` \tab probability mass left out from tails if integration needs to be restricted in range \cr
#' }
#'
#' @section Version History:
#'
#' See `NEWS.md` file.
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
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
#' @importFrom utils capture.output modifyList
## usethis namespace: end
"_PACKAGE"
