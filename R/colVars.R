#' Fast column-wise calculation of unbiased variances
#' @keywords internal
colVars <- function(a) {
    n <- dim(a)[[1]]
    c <- dim(a)[[2]]
    return(.colMeans(((a - matrix(.colMeans(a, n, c), nrow = n, ncol = c, byrow = TRUE)) ^ 2), n, c) * n / (n - 1))
}

