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
#' The main function of the package is \code{\link{gMAP}}. See it's
#' help page for a detailed description of the statistical model.
#'
#'
#' @section Global Options:
#'
#' \tabular{lcl}{
#' Option \tab Default \tab Description \cr
#' \code{RBesT.MC.warmup} \tab 2000 \tab MCMC warmup iterations \cr
#' \code{RBesT.MC.iter} \tab 6000 \tab total MCMC iterations \cr
#' \code{RBesT.MC.chains} \tab 4 \tab MCMC chains\cr
#' \code{RBesT.MC.thin} \tab 4 \tab MCMC thinning \cr
#' \code{RBesT.MC.control} \tab \code{list(adapt_delta=0.99,} \tab sets \code{control} argument for Stan call\cr
#'  \tab \code{stepsize=0.01,} \tab \cr
#'  \tab \code{max_treedepth=20)} \tab \cr
#' \code{RBesT.MC.ncp} \tab 1 \tab parametrization: 0=CP, 1=NCP, 2=Automatic  \cr
#' \code{RBesT.MC.init} \tab 1 \tab range of initial uniform [-1,1] is the default  \cr
#' \code{RBesT.MC.rescale} \tab \code{TRUE} \tab Automatic rescaling of raw parameters  \cr
#' \code{RBesT.verbose} \tab \code{FALSE} \tab requests outputs to be more verbose\cr
#' \code{RBesT.integrate_args} \tab \code{list(lower=-Inf,} \tab arguments passed to \code{integrate} for\cr
#'  \tab \code{upper=Inf,} \tab intergation of densities\cr
#' \tab \code{rel.tol=.Machine$double.eps^0.25,} \tab \cr
#' \tab \code{abs.tol=.Machine$double.eps^0.25,} \tab \cr
#' \tab \code{subdivisions=1E3)} \tab \cr
#' }
#'
#' @section Version History:
#'
#' See \code{NEWS} file.
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#'
#' @name RBesT-package
#' @aliases RBesT
#' @docType package
#' @useDynLib RBesT, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom RcppParallel RcppParallelLibs CxxFlags
#' @importFrom rstan sampling extract get_sampler_params summary
#' @import rstantools
#' @import stats
#' @importFrom utils capture.output modifyList
#' @importFrom matrixStats rowLogSumExps colLogSumExps colSums2 rowMins rowRanks logSumExp
#' @import assertthat
#' @import mvtnorm
#' @import ggplot2
#' @import Formula
#' @import checkmate
#' @import abind
NULL
