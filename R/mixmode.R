mixmode <- function(mix,  interval=0.99, verbose=FALSE) {
    tol <- .Machine$double.eps^0.25
    digits <- floor(abs(log10(tol)))
    compmode <- function(comp) {
        mixComp <- mix[[comp]]
        mixComp[1,] <- 1
        qlow <- (1-interval)/2
        qup <- 1-qlow
        optimise(function(x) dmix(mix, x, log=TRUE), qmix(mixComp, c(qlow, qup)), maximum = TRUE, tol=tol)$maximum
    }
    ## find around each component the maximum as otherwise optimise
    ## can run into a local extremum
    extrema <- sapply(seq(ncol(mix)), compmode)
    res <- extrema[which.max(dmix(mix, extrema))]
    ## identify distinct modes
    ind <- duplicated(signif(extrema, digits))
    attr(res, "modes") <- extrema[!ind]
    if(verbose) {
        cat("Locations:",extrema[!ind], "\n")
        cat("Density  :",dmix(mix, extrema[!ind]), "\n")
    }
    res
}
