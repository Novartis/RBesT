#' Extract log likelihood from fitted EM objects
#'
#' @keywords internal
#' 
#' @export
logLik.EM <- function(object, ...) {
    val <- attr(object, "lli")
    attr(val, "df") <-  attr(object, "df")
    attr(val, "nobs") <-  attr(object, "nobs")
    class(val) <- "logLik"
    val
}
