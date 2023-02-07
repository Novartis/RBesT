#' Summarize Arrays
#'
#' The function calculates summary statistics from arbitrary arrays.
#'
#' @param x Object to summarize which can be a numerical vector, matrix or a multi-dimensional array
#' @param min.max Enables to include minimum and maximum in the output.
#' @param n.sim Enables to include the number of observations in the output.
#' @param probs Quantiles to output.
#' @param margin Margin of the input array over which the summary function is applied.
#'
#' @details The function calculates by default the mean, standard
#' deviation and the specified qantiles which are by default the
#' median and the 95% interval.
#'
#' If a mulit-dimensional array is specified as \code{x}, then the
#' function will by default calculate the summaries over the margin of
#' the largest dimension. For the case of a vector and a matrix, the
#' function will transpose the results for better readabiliy.
#'
#' @keywords internal
#'
`SimSum` <-
function( x, min.max=FALSE, n.sim=FALSE, probs=c(0.025,0.5,0.975), margin=ifelse(is.null(dim(x) | length(dim(x)) == 1), 2, length(dim(x))) )
{
    #	Version 1.1, 16-Oct-2014
    if(is.null(dim(x)) | length(dim(x)) == 1)
       x <- matrix(x, ncol=1)

    sfun <- function(r) {
        ex <- if(min.max) c(min=min(r), max=max(r)) else NULL
        N  <- if(n.sim)   c(nsim=length(r))         else NULL
        c( mean=mean(r), sd=sd(r), quantile(r,probs), ex, N )
    }

    sim.out <- apply(x, margin, sfun)

    ## to ensure compatibility with old versions of the function, a
    ## transpose is needed for the standard case of margin being 2
    if(is.matrix(x) && margin == 2)
        sim.out <- t(sim.out)

    return(sim.out)
}

