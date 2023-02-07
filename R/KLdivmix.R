## Kullback-Leibler distance between mixtures
KLdivmix <- function(mixRef, mixTest) {
    interval <- support(mixRef)
    if(!all(interval == support(mixRef)))
        warning("Support of mixRef and mixTest do not match.")
    ## note: setting stop.on.error to FALSE manages to avoid boundary
    ## value issues
    integrate(function(x) { dmix(mixRef, x) * (dmix(mixRef,x,log=TRUE) - dmix(mixTest,x,log=TRUE)) }, interval[1], interval[2], stop.on.error=FALSE )$value
}
