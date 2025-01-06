#' @name mixmvnorm
#'
#' @title Multivariate Normal Mixture Density
#'
#' @description The multivariate normal mixture density and auxiliary
#'     functions.
#'
#' @param ... List of mixture components.
#' @param sigma Reference covariance.
#' @param param Determines how the parameters in the list are
#'     interpreted. See details.
#' @param object Multivariate normal mixture object.
#'
#' @details Each entry in the \code{...} argument list is a numeric
#'     vector defining one component of the mixture multivariate
#'     normal distribution. The first entry of the component defining
#'     vector is the weight of the mixture component followed by the
#'     vector of means in each dimension and finally a specification
#'     of the covariance matrix, which depends on the chosen
#'     parametrization. The covariance matrix is expected to be given
#'     as numeric vector in a column-major format, which is standard
#'     conversion applied to matrices by the vector concatenation
#'     function \code{\link[base:c]{c}}. Please refer to the examples
#'     section below.
#'
#' Each component defining vector can be specified in different ways
#' as determined by the \code{param} option:
#'
#' \describe{
#' \item{ms}{Mean vector and covariance matrix \code{s}. Default.}
#' \item{mn}{Mean vector and number of observations. \code{n} determines the covariance for each component via the relation \eqn{\Sigma/n} with \eqn{\Sigma} being the known reference covariance.}
#' }
#'
#' The reference covariance \eqn{\Sigma} is the known covariance in
#' the normal-normal model (observation covariance). The function
#' \code{sigma} can be used to query the reference covariance and may
#' also be used to assign a new reference covariance, see examples
#' below. In case \code{sigma} is not specified, the user has to
#' supply \code{sigma} as argument to functions which require a
#' reference covariance.
#'
#' @family mixdist
#'
#' @return Returns a multivariate normal mixture with the specified
#'     mixture components.
#'
#' @examples
#'
#' S <- diag(c(1, 2)) %*% matrix(c(1, 0.5, 0.5, 1), 2, 2) %*% diag(c(1, 2))
#' mvnm1 <- mixmvnorm(
#'   rob = c(0.2, c(0, 0), diag(c(5, 5))),
#'   inf = c(0.8, c(0.5, 1), S / 10), sigma = S
#' )
#'
#' print(mvnm1)
#' summary(mvnm1)
#'
#' set.seed(657846)
#' mixSamp1 <- rmix(mvnm1, 500)
#' colMeans(mixSamp1)
#'
NULL

#' @rdname mixmvnorm
#' @export
mixmvnorm <- function(..., sigma, param = c("ms", "mn")) {
  ## length of first mean vector determines dimension
  mix <- mixdist3(...)
  dim_labels <- rownames(mix)
  param <- match.arg(param)
  Nc <- ncol(mix)
  n <- colnames(mix)
  if (param == "ms") {
    ## mean vector & covariance parametrization
    l <- nrow(mix)
    p <- (sqrt(1 + 4 * (l - 1)) - 1) / 2
    assert_integerish(p, lower = 1, any.missing = FALSE, len = 1)
    ## in this case we expect c(weight, mean, as.numeric(cov))
    mix <- do.call(
      mixdist3,
      lapply(
        1:Nc,
        function(co) c(mix[1, co], mvnorm(mix[2:(p + 1), co], matrix(mix[(1 + p + 1):(1 + p + p^2), co], p, p)))
      )
    )
  }
  if (param == "mn") {
    ## mean vector & number of observations
    l <- nrow(mix)
    p <- l - 2
    assert_integerish(p, lower = 1, any.missing = FALSE, len = 1)
    assert_matrix(sigma, any.missing = FALSE, nrows = p, ncols = p)
    mix <- do.call(mixdist3, lapply(
      1:Nc,
      function(co) {
        assert_numeric(mix[l, co], lower = 0, finite = TRUE, any.missing = FALSE)
        c(mix[1, co], mvnorm(mix[2:(p + 1), co], sigma / mix[l, co]))
      }
    ))
  }
  colnames(mix) <- n
  p <- mvnormdim(mix[-1, 1])
  if (is.null(dim_labels)) {
    dim_labels <- as.character(1:p)
  } else {
    dim_labels <- dim_labels[2:(p + 1)]
  }
  rownames(mix) <- c("w", mvnorm_label(mix[-1, 1], dim_labels))
  if (!missing(sigma)) {
    assert_matrix(sigma, any.missing = FALSE, nrows = p, ncols = p)
    colnames(sigma) <- rownames(sigma) <- dim_labels
    attr(mix, "sigma") <- sigma
  }
  class(mix) <- c("mvnormMix", "mix")
  likelihood(mix) <- "mvnormal"
  mix
}

#' @keywords internal
mvnorm <- function(mean, sigma) {
  ## TODO: far more checks!!! Allow to pass in directly a cholesky factor
  assert_numeric(mean, finite = TRUE, any.missing = FALSE)
  p <- length(mean)
  assert_matrix(sigma, any.missing = FALSE, nrows = p, ncols = p)
  rho <- cov2cor(sigma)
  s <- sqrt(diag(sigma))
  mvn <- c(mean, s, rho[lower.tri(rho)])
  names(mvn) <- mvnorm_label(mvn)
  mvn
}

#' @keywords internal
mvnormdim <- function(mvn) {
  p <- (-3 + sqrt(9 + 8 * length(mvn))) / 2
  assert_integerish(p, lower = 1, any.missing = FALSE, len = 1)
  as.integer(p)
}

#' @keywords internal
mvnorm_dim_labels <- function(mvn) {
  p <- mvnormdim(mvn)
  n <- names(mvn)
  if (is.null(n)) {
    return(as.character(1:p))
  }
  return(gsub("\\]$", "", gsub("^[^\\[]*\\[", "", n[1:p])))
}

#' @keywords internal
mvnorm_label <- function(mvn, dim_labels) {
  p <- mvnormdim(mvn)
  if (missing(dim_labels)) {
    dim_labels <- mvnorm_dim_labels(mvn)
  }
  if (p > 1) {
    Rho_labs_idx <- outer(dim_labels, dim_labels, paste, sep = ",")[lower.tri(diag(p))]
    lab <- c(paste0("m[", dim_labels, "]"), paste0("s[", dim_labels, "]"), paste0("rho[", Rho_labs_idx, "]"))
  } else {
    lab <- c(paste0("m[", dim_labels, "]"), paste0("s[", dim_labels, "]"))
  }
  lab
}

#' @keywords internal
mvnormsigma <- function(mvn) {
  p <- mvnormdim(mvn)
  n <- length(mvn)
  s <- mvn[(p + 1):(2 * p)]
  Rho <- diag(p)
  Rho[lower.tri(Rho)] <- mvn[(2 * p + 1):n]
  Rho[upper.tri(Rho)] <- t(Rho)[upper.tri(Rho)]
  diag(s, nrow = p) %*% Rho %*% diag(s, nrow = p)
}


#' @rdname mixmvnorm
#' @method print mvnormMix
#' @param x The mixture to print
#' @export
print.mvnormMix <- function(x, ...) {
  cat("Multivariate normal mixture\n")
  cat(paste0("Outcome dimension: ", mvnormdim(x[-1, 1]), "\n"))
  if (!is.null(sigma(x))) {
    cat("Reference covariance:\n")
    print(sigma(x), ...)
  }
  NextMethod()
}

#' @rdname mixmvnorm
#' @method summary mvnormMix
#' @export
summary.mvnormMix <- function(object, ...) {
  w <- object[1, ]
  p <- mvnormdim(object[-1, 1])
  m <- object[2:(p + 1), , drop = FALSE]
  Nc <- ncol(object)
  mmix <- rowSums(sweep(m, 2, w, "*"))
  ## Cov(x,x) = E(x x') - E(x) E(x')
  ## E(x x') = Sigma + m m' (see matrix cookbook eq 377)
  if (Nc == 1) {
    S <- w[1] * mvnormsigma(object[-1, 1])
  } else {
    S <- -1 * tcrossprod(mmix)
    for (i in 1:Nc) {
      S <- S + w[i] * (mvnormsigma(object[-1, i]) + tcrossprod(unname(m[, i, drop = FALSE])))
    }
  }
  rownames(S) <- colnames(S) <- names(mmix) <- mvnorm_dim_labels(object[-1, 1])
  list(mean = mmix, cov = S)
}

#' @rdname mixmvnorm
#' @method sigma mvnormMix
#' @export
#' @export sigma
sigma.mvnormMix <- function(object, ...) {
  attr(object, "sigma")
}

#' @keywords internal
is_mixmv <- function(mix) {
  inherits(mix, "mvnormMix")
}

#' @keywords internal
mvnorm_extract_dim <- function(mix, sub) {
  Nc <- ncol(mix)
  sub_comp <- list()
  p <- mvnormdim(mix[-1, 1])
  assert_numeric(sub, lower = 1, upper = p, any.missing = FALSE)
  for (i in seq_len(Nc)) {
    sub_comp[[i]] <- c(mix["w", i], mix[1 + sub, i], mvnormsigma(mix[-1, i])[sub, sub, drop = FALSE])
  }
  if (!is.null(sigma(mix))) {
    sub_comp$sigma <- sigma(mix)
  }
  do.call(mixmvnorm, sub_comp)
}
