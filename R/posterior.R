# Code adapted from
# https://raw.githubusercontent.com/paul-buerkner/brms/refs/heads/master/R/posterior.R

#'
#' Transform `gMAP` to `draws` objects
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Transform a `gMAP` object to a format supported by the
#' \pkg{posterior} package.
#'
#' @aliases as_draws as_draws_matrix as_draws_array as_draws_df
#' @aliases as_draws_rvars as_draws_list
#'
#' @param x A `gMAP` object.
#' @param variable A character vector providing the variables to
#'   extract.  By default, all variables are extracted.
#' @param regex Logical; Should variable be treated as a (vector of)
#'   regular expressions? Any variable in `x` matching at least
#'   one of the regular expressions will be selected. Defaults to
#'   `FALSE`.
#' @param inc_warmup Should warmup draws be included? Defaults to
#'   `FALSE`.
#' @param ... Arguments passed to individual methods (if applicable).
#'
#' @details To subset iterations, chains, or draws, use the
#'   [posterior::subset_draws()] method after
#'   transforming the input object to a `draws` object.
#'
#' The function is experimental as the set of exported posterior
#' variables are subject to updates.
#'
#' @seealso [posterior::draws()]
#'   [posterior::subset_draws()]
#'
#' @template example-start
#' @examples
#'
#' set.seed(34563)
#' map_AS <- gMAP(cbind(r, n - r) ~ 1 | study,
#'   family = binomial,
#'   data = AS,
#'   tau.dist = "HalfNormal", tau.prior = 1,
#'   beta.prior = 2
#' )
#'
#' post_AS <- as_draws(map_AS)
#'
#' @template example-stop
#'
#' @name draws-RBesT
NULL

#' @rdname draws-RBesT
#' @importFrom posterior as_draws
#' @method as_draws gMAP
#' @export
#' @export as_draws
as_draws.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  # draws_list is the fastest format to convert to at the moment
  .as_draws_conversion(
    x$fit,
    as_draws_list,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-RBesT
#' @importFrom posterior as_draws_matrix
#' @method as_draws_matrix gMAP
#' @export
#' @export as_draws_matrix
as_draws_matrix.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x$fit,
    as_draws_matrix,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-RBesT
#' @importFrom posterior as_draws_array
#' @method as_draws_array gMAP
#' @export
#' @export as_draws_array
as_draws_array.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x$fit,
    as_draws_array,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-RBesT
#' @importFrom posterior as_draws_df
#' @method as_draws_df gMAP
#' @export
#' @export as_draws_df
as_draws_df.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x$fit,
    as_draws_df,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-RBesT
#' @importFrom posterior as_draws_list
#' @method as_draws_list gMAP
#' @export
#' @export as_draws_list
as_draws_list.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x$fit,
    as_draws_list,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-RBesT
#' @importFrom posterior as_draws_rvars
#' @method as_draws_rvars gMAP
#' @export
#' @export as_draws_rvars
as_draws_rvars.gMAP <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x$fit,
    as_draws_rvars,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

# in stanfit objects draws are stored in a draws_list-like format
# so converting from there will be most efficient
# may be removed once rstan supports posterior natively
#' @keywords internal
.as_draws_conversion <- function(
  x,
  draws_converter,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  stopifnot(.is.stanfit(x))
  inc_warmup <- .as_one_logical(inc_warmup)
  if (!length(x@sim$samples)) {
    .stop2("The model does not contain posterior draws.")
  }
  ## since rstan::extract returns an array we tell posterior to go via
  ## arrays directly
  out <- draws_converter(as_draws_array(rstan::extract(
    x,
    permuted = FALSE,
    inc_warmup = inc_warmup
  )))
  # subset variables
  subset_draws(out, variable = variable, regex = regex)
}

#' @keywords internal
.is.stanfit <- function(x) {
  inherits(x, "stanfit")
}

# coerce 'x' to a single logical value
#' @keywords internal
.as_one_logical <- function(x, allow_na = FALSE) {
  s <- substitute(x)
  x <- as.logical(x)
  if (length(x) != 1L || anyNA(x) && !allow_na) {
    s <- .deparse0(s, max_char = 100L)
    .stop2("Cannot coerce '", s, "' to a single logical value.")
  }
  x
}

#' @keywords internal
.stop2 <- function(...) {
  stop(..., call. = FALSE)
}

#' @keywords internal
collapse <- function(..., sep = "") {
  paste(..., sep = sep, collapse = "")
}

# combine deparse lines into one string
# since R 4.0 we also have base::deparse1 for this purpose
#' @keywords internal
.deparse0 <- function(x, max_char = NULL, ...) {
  out <- collapse(deparse(x, ...))
  if (isTRUE(max_char > 0)) {
    out <- substr(out, 1L, max_char)
  }
  out
}
