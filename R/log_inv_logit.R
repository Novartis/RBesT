#' Numerically stable log of the inv_logit function
#' @keywords internal
log_inv_logit <- function(mat) {
    ##- ifelse(is.finite(mat) & (mat < 0), log1p(exp(mat)) - mat, log1p(exp(-mat)))
    ##idx <- is.finite(mat) & (mat < 0)
    idx <- mat < 0
    mat[idx] <- mat[idx] - log1p(exp(mat[idx]))
    mat[!idx] <- -1*log1p(exp(-mat[!idx]))
    mat
}
