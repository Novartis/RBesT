#' Scrambles the order of a mcmc array object for usage as a mcmc
#' sample. It is advisable to set order once per mcmc run, otherwise
#' correlations in the mcmc sample will be lost.
#' @keywords internal
chains2sample <- function(chains, order, drop=TRUE) {
    d <- dim(chains)
    N <- prod(d[2:3])
    if(missing(order))
        order <- sample.int(N)
    matrix(chains, nrow=d[1], ncol=N)[,order,drop=drop]
}
