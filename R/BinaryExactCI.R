#' Exact Confidence interval for Binary Proportion
#'
#' This function calculates the exact confidendence interval for a
#' response rate presented by \eqn{n} and \eqn{r}.
#'
#' @param r Number of success or responder 
#' @param n Sample size
#' @param alpha  confidence level
#' @param drop Determines if \code{\link{drop}} will be called on the result
#' 
#' @details
#' Confidence intervals are obtained by a procedure first given in 
#' Clopper and Pearson (1934). This guarantees that the confidence 
#' level is at least (1-\eqn{\alpha}).
#'
#' Details can be found in the publication listed below.
#' 
#' @return 100 (1-\eqn{\alpha})\% exact confidence interval for given
#' response rate
#'
#' @references Clopper, C. J. & Pearson, E. S. The use of confidence or
#' fiducial limits illustrated in the case of the binomial. Biometrika 1934. 
#'  
#' @examples
#' BinaryExactCI(3,20,0.05)
#' 
#' @export
BinaryExactCI <- function(r, n, alpha=0.05, drop=TRUE) {
    alpha2 <- alpha/2
    Low <- alpha2
    High <- 1-alpha2
    
    pLow <- qbeta( Low, r+(r==0), n-r+1)
    pHigh <- qbeta( High, r+1, n-r+((n-r)==0))

    nms <- c( paste(round(100*Low,1),"%",sep=""),paste(round(100*High,1),"%",sep="") )

    CI <- cbind(pLow,pHigh)
    colnames(CI) <- nms

    if(drop) CI <- drop(CI)
    
    return( CI )
}


