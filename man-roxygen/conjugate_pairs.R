#' @section Supported Conjugate Prior-Likelihood Pairs:
#' 
#' \tabular{lccc}{
#' \strong{Prior/Posterior} \tab \strong{Likelihood} \tab \strong{Predictive} 
#'  \tab \strong{Summaries} \cr
#' Beta \tab Binomial \tab Beta-Binomial \tab \code{n}, \code{r} \cr
#' Normal \tab Normal (\emph{fixed \eqn{\sigma}}) \tab Normal \tab \code{n}, \code{m}, \code{se}  \cr
#' Gamma \tab Poisson \tab Gamma-Poisson \tab  \code{n}, \code{m} \cr
#' Gamma \tab Exponential \tab Gamma-Exp (\emph{not supported}) \tab \code{n}, \code{m}
#' }
#' 
