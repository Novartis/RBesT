#' Transform Densities with a link function
#'
#' One-to-one transforms (mixture) of densities using a link function.
#'
#' @param object Mixture density to apply link to.
#' @param value Link.
#'
#' Note: link functions are assumed to be order preserving, i.e. if
#' x_1 < x_2 holds, then link(x_1) < link(x_2).
#'
#' @keywords internal
'dlink<-' <- function(object, value) {
    if(is.dlink(value))
        trans <- value
    else {
        trans <- match.fun(value)()
    }
    assert_that(is.dlink(trans))
    attr(object, "link") <- trans
    object
}

dlink <- function(object) {
    attr(object, "link")
}

dlink_new <- function(name, link, inv, Jinv_orig, lJinv_orig, lJinv_link) {
    if (is.character(link))
        link <- match.fun(link)
    if (is.character(inv))
        inv <- match.fun(inv)
    if (is.character(Jinv_orig))
        Jinv_orig <- match.fun(Jinv_orig)
    if (is.character(lJinv_orig))
        lJinv_orig <- match.fun(lJinv_orig)
    if (is.character(lJinv_link))
        lJinv_link <- match.fun(lJinv_link)

    structure(list(name=name, link=link, invlink=inv,
                   Jinv_orig=Jinv_orig,
                   lJinv_orig=lJinv_orig,
                   lJinv_link=lJinv_link), class="dlink")
}

identity_dlink <- dlink_new("identity",
                            identity, identity,
                            Curry(fill, value=1), Curry(fill, value=0), Curry(fill, value=0))

logit_Jinverse_orig <- function(mu) mu * (1-mu)
logit_lJinverse_orig <- function(mu) log(mu) + log1p(-mu)
logit_lJinverse_link <- function(l) log_inv_logit(l) + log_inv_logit(-l)

logit_dlink <- dlink_new("logit",
                         binomial()$linkfun, binomial()$linkinv,
                         logit_Jinverse_orig, logit_lJinverse_orig, logit_lJinverse_link)

log_dlink <- dlink_new("log",
                       log, exp,
                       identity, log, identity)

link_map <- list(identity=identity_dlink,
                 logit=logit_dlink,
                 log=log_dlink)

canonical_dlink <- function(mix) {
    if(inherits(mix, "betaMix")) {
        dlink(mix) <- logit_dlink
        return(mix)
    }
    if(inherits(mix, "gammaMix")) {
        dlink(mix) <- log_dlink
        return(mix)
    }
    return(mix)
}

#' Fill numeric objects
#'
#' Returns the numeric input object with the value given and respects
#' dimensionalty and type of input.
#'
#' @param x Input numeric object.
#' @param value Value filled.
#'
#' @keywords internal
fill <- function(x, value) {
    cl <- class(x)
    ax <- array(value, dim(as.array(x)))
    class(ax) <- class(x)
    attributes(ax) <- attributes(x)
    ax
}

is.dlink <- function(x)
    inherits(x, "dlink")

is.dlink_identity <- function(x)
    is.dlink(x) & x$name == "identity"

print.dlink <- function(x, ...)
    cat("Link:", x$name, "\n")
