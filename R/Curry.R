#' Functional programming utilities
#' 
#' function from functional
#' 
#' @keywords internal
Curry <- function (FUN, ...) 
{
    .orig = list(...)
    function(...) do.call(FUN, c(.orig, list(...)))
}
