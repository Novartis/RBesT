#' Find root of univariate function of integers
#'
#' Uses a bisectioning algorithm to search the give interval for a
#' change of sign and returns the integer which is closest to 0.
#'
#' @keywords internal
uniroot_int <- function(f, interval, ...,
                        f.lower=f(interval[1], ...),
                        f.upper=f(interval[2], ...),
                        maxIter=1000) {
    lo <- interval[1]
    hi <- interval[2]

    assert_that(interval[1] < interval[2])

    fleft <- f.lower
    fright <- f.upper

    if(f.lower*f.upper > 0) {
        ##warning("Minimum not in range!")
        ##return(ifelse(abs(flo) < abs(fhi), lo, hi))
        return(numeric())
    }

    iter <- 0
    while ((hi-lo) > 1 & iter < maxIter) {
        mid <- floor((lo + hi) / 2)
        fmid <- f(mid, ...)
        if (f.lower * fmid < 0) { hi <- mid; fright <- fmid; }
        else if (f.upper * fmid < 0) { lo <- mid; fleft <- fmid; }
        iter <- iter + 1
    }
    if(iter == maxIter)
        warning("Maximum number of iterations reached.")
    return(ifelse(abs(fleft) < abs(fright), lo, hi))
}

uniroot_int.all <- function (f, interval, maxIter=1000, n = 100, ...) 
{
    assert_that(interval[1] < interval[2])

    xseq <- round(seq(interval[1], interval[2], len = n + 1))
    xseq <- xseq[!duplicated(xseq)]
    nu <- length(xseq) - 1
    mod <- f(xseq, ...)
    Equi <- xseq[which(mod == 0)]
    ss <- mod[1:nu] * mod[2:(nu + 1)]
    print(ss)
    ii <- which(ss < 0)
    print(ii)
    print(xseq[c(ii, ii[length(ii)] + 1)])
    for (i in ii) Equi <- c(Equi, uniroot_int(f, c(xseq[i], xseq[i + 1]), ...,
                                              maxIter=maxIter))
    return(Equi)
}
