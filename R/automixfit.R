#' Automatic Fitting of Mixtures of Conjugate Distributions to a Sample
#'
#' Fitting a series of mixtures of conjugate distributions to a
#' \code{sample}, using Expectation-Maximization (EM). The number of
#' mixture components is specified by the vector \code{Nc}. First a
#' \code{Nc[1]} component mixture is fitted, then a \code{Nc[2]}
#' component mixture, and so on. The mixture providing the best AIC
#' value is then selected.
#'
#' @param sample Sample to be fitted by a mixture distribution.
#' @param Nc Vector of mixture components to try out (default \code{seq(1,4)}).
#' @param k Penalty parameter for AIC calculation (default 6)
#' @param thresh The procedure stops if the difference of subsequent AIC values
#' is smaller than this threshold (default -Inf). Setting the threshold to 0
#' stops \code{automixfit} once the AIC becomes worse.
#' @param verbose Enable verbose logging.
#' @param ... Further arguments passed to \code{\link{mixfit}},
#' including \code{type}.
#'
#' @details The \code{type} argument specifies the distribution of
#' the mixture components, and can be a normal, beta or gamma
#' distribution.
#'
#' The penalty parameter \code{k} is 2 for the standard AIC
#' definition. \emph{Collet (2003)} suggested to use values in the
#' range from 2 to 6, where larger values of \code{k} penalize more
#' complex models. To favor mixtures with fewer components a value of
#' 6 is used as default.
#'
#' @return As result the best fitting mixture model is returned,
#' i.e. the model with lowest AIC. All other models are saved in the
#' attribute \code{models}.
#'
#' @references
#' Collet D.
#' \emph{Modeling Survival Data in Medical Research}.
#' 2003; Chapman and Hall/CRC.
#'
#' @examples
#' # random sample of size 1000 from a mixture of 2 beta components
#' bm <- mixbeta(beta1=c(0.4, 20, 90), beta2=c(0.6, 35, 65))
#' bmSamp <- rmix(bm, 1000)
#'
#' # fit with EM mixture models with up to 10 components and stop if
#' # AIC increases
#' bmFit <- automixfit(bmSamp, Nc=1:10, thresh=0, type="beta")
#' bmFit
#'
#' # advanced usage: find out about all discarded models
#' bmFitAll <- attr(bmFit, "models")
#'
#' sapply(bmFitAll, AIC, k=6)
#'
#'
#' @export
automixfit <- function(sample, Nc=seq(1, 4), k=6, thresh=-Inf, verbose=FALSE, ...) {
    if("gMAPpred" %in% class(sample)) {
        stop("Not yet supported.")
    }
    assert_that(all(diff(Nc) >=1 ))
    models <- list()
    ic <- Inf
    aic <- c()
    for(i in seq_along(Nc)) {
        curNc <- Nc[i]
        icLast <- ic
        if(!verbose)
            suppressMessages(run <- mixfit(sample, Nc=curNc, verbose=verbose, ...))
        else
            run <- mixfit(sample, Nc=curNc, verbose=verbose, ...)
        ic <- AIC(run, k=k)
        aic <- c(aic, ic)
        delta <- icLast - ic
        if(verbose)
            message("Components", curNc, ": AIC", ic, "; deltaAIC =", delta, "\n")
        models <- c(models, list(run))
        if(delta < thresh)
            break
    }
    names(models) <- Nc[1:i]
    o <- order(aic)
    models <- models[o]
    models
    bestfit <- models[[1]]
    attr(bestfit, "models") <- models
    bestfit
}

