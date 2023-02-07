#' Plot mixture distributions
#'
#' @param x mixture distribution
#' @param prob  defining lower and upper percentile of x-axis. Defaults to the 99\% central probability mass.
#' @param fun function to plot which can be any of \code{dmix}, \code{qmix} or \code{pmix}.
#' @param log log argument passed to the function specified in \code{fun}.
#' @param comp for the density function this can be set to \code{TRUE}
#' which will display colour-coded each mixture component of the
#' density in addition to the density.
#' @param size controls the linesize in plots.
#' @param ... extra arguments passed on to the \code{\link[ggplot2]{qplot}} call.
#'
#' @details Plot function for mixture distribution objects. It shows
#' the density/quantile/cumulative distribution (corresponds to
#' \code{d/q/pmix} function) for some specific central probability
#' mass defined by \code{prob}. By default the x-axis is chosen to
#' show 99\% of the probability density mass.
#'
#' @template plot-help
#'
#' @return
#' A \code{\link[ggplot2]{ggplot}} object is returned.
#'
#' @family mixdist
#' @examples
#' # beta with two informative components
#' bm <- mixbeta(inf=c(0.5, 10, 100), inf2=c(0.5, 30, 80))
#' plot(bm)
#' plot(bm, fun=pmix)
#'
#' # for customizations of the plot we need to load ggplot2 first
#' library(ggplot2)
#'
#' # show a histogram along with the density
#' plot(bm) + geom_histogram(data=data.frame(x=rmix(bm, 1000)),
#'                           aes(y=..density..), bins=50, alpha=0.4)
#'
#' \donttest{
#' # note: we can also use bayesplot for histogram plots with a density ...
#' library(bayesplot)
#' mh <- mcmc_hist(data.frame(x=rmix(bm, 1000)), freq=FALSE) +
#'          overlay_function(fun=dmix, args=list(mix=bm))
#' # ...and even add each component
#' for(k in 1:ncol(bm))
#'   mh <- mh + overlay_function(fun=dmix, args=list(mix=bm[[k]]), linetype=I(2))
#' print(mh)
#' }
#'
#' # normal mixture
#' nm <- mixnorm(rob=c(0.2, 0, 2), inf=c(0.8, 6, 2), sigma=5)
#' plot(nm)
#' plot(nm, fun=qmix)
#'
#' # obtain ggplot2 object and change title
#' pl <- plot(nm)
#' pl + ggtitle("Normal 2-Component Mixture")
#'
#' @method plot mix
#' @export
plot.mix <- function(x, prob=0.99, fun=dmix, log=FALSE, comp=TRUE, size=1.25, ...) {
    funStr <- deparse(substitute(fun))
    if(length(prob) == 1) {
        plow <- (1-prob)/2
        pup <- 1-plow
        if(funStr != "qmix") {
            interval <- qmix(x, c(plow, pup))
        } else {
            interval <- c(plow, pup)
        }
    } else {
        plow <- prob[1]
        pup <- prob[2]
        interval <- prob
    }
    assert_that(plow < pup)
    assert_that(interval[1] < interval[2])
    fun <- match.fun(fun)
    discrete <- ifelse(all(is.integer(interval)), TRUE, FALSE )
    if(discrete) {
        plot_fun <- function(x, ...) fun(floor(x), ...)
        plot_geom <- "step"
    } else {
        plot_fun <- fun
        plot_geom <- "line"
    }
    n_fun <- 501

    num_comp <- ncol(x)
    pl <- ggplot(data.frame(x=interval), aes(x=x)) +
        stat_function(geom=plot_geom, fun = plot_fun, args=list(mix=x, log=log), n=n_fun, size=size) +
        bayesplot::bayesplot_theme_get()

    if(funStr=="dmix") {
        pl <- pl + ylab("density") + xlab("parameter")
    } else if(funStr=="pmix") {
        pl <- pl + ylab("cumulative density") + xlab("quantile")
    } else if(funStr=="qmix") {
        pl <- pl + ylab("quantile") + xlab("cumulative density")
    }
    if(funStr=="dmix" & comp) {
        comp_df <- list()
        for(i in seq(num_comp)) {
            pl <- pl + stat_function(geom=plot_geom, mapping=aes_(colour=factor(i)), fun=plot_fun, args=list(mix=x[[i]], log=log), n=n_fun, linetype=I(2), size=size)
        }
        pl <- pl + scale_colour_manual("Comp. [%]", values=2:(num_comp+1), labels=paste0(colnames(x), " ", format(100*x[1,],digits=1,nsmall=1)))
    }
    pl
}
