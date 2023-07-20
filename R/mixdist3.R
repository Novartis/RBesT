#' Utility function to instantiate 2 parameter mixture densities.
#'
#' @keywords internal
mixdist3 <- function(...) {
    args <- list(...)
    Nc <- length(args)
    l <- sapply(args, length)
    if(!all(l >= 3) | !all(l == l[1]))
        stop("All components must have equal number of parameters.")
    res <- do.call(cbind, args)
    if(is.null(names(args)))
        colnames(res) <- paste("comp", seq(Nc), sep="")
    norm <- sum(res[1,])
    if(norm != 1) {
        ## only issue a warning if difference appears to be a real
        ## user error, otherwise just silently rescale since we are
        ## correcting floating point arithmetic unless RBesT is asked
        ## to be verbose
        if(getOption("RBesT.verbose", FALSE) | abs(norm - 1) > 1E-4)
            warning("Weights do not sum to 1. Rescaling accordingly.")
        res[1,] <- res[1,]/norm
    }
    ## assign the default identity transform
    dlink(res) <- identity_dlink
    res
}
