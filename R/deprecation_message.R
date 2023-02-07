deprecated <- function(what, replaced) {
    message(what, " is deprecated and will be removed in a future release.")
    if(!missing(replaced))
        message("Please use instead ", replaced, ".")
}
