#' Forest Plot
#'
#' Creates a forest plot for \code{\link{gMAP}} analysis objects.
#'
#' @param x \code{\link{gMAP}} object.
#' @param prob confidence interval width and probability mass of credible intervals.
#' @param est can be set to one of \code{both} (default), \code{MAP}, \code{Mean} or \code{none}. Controls which model estimates are to be included.
#' @param model controls which estimates are displayed per study. Either \code{stratified} (default), \code{both} or \code{meta}.
#' @param point_est shown point estimate. Either \code{median} (default) or \code{mean}.
#' @param size controls point and linesize.
#' @param alpha transparency of reference line. Setting \code{alpha=0}
#' suppresses the reference line.
#'
#' @details The function creates a forest plot suitable for
#' \code{\link{gMAP}} analyses. Note that the Meta-Analytic-Predictive
#' prior is included by default in the plot as opposed to only showing
#' the estimated model mean. See the examples below to obtain standard
#' forest plots.
#'
#' Also note that the plot internally flips the x and
#' y-axis. Therefore, if you want to manipulate the x-axis, you have
#' to give commands affecting the y-axis (see examples).
#'
#' @template plot-help
#'
#' @return The function returns a \pkg{ggplot2} plot object.
#'
#' @seealso \code{\link{gMAP}}
#'
#' @examples
#' # we consider the example AS MAP analysis
#' example(AS)
#'
#' # default forest plot for a gMAP analysis
#' forest_plot(map_AS)
#'
#' # standard forest plot (only stratified estimate and Mean)
#' forest_plot(map_AS, est=c("Mean"), model="stratified")
#'
#' # to further customize these plots, first load bayesplot and ggplot2
#' library(bayesplot)
#' library(ggplot2)
#'
#' # to make plots with red colors, big fonts for presentations, suppress
#' # the x axis label and add another title (with a subtitle)
#' color_scheme_set("red")
#' theme_set(theme_default(base_size=16))
#' forest_plot(map_AS, size=2) +
#'    yaxis_title(FALSE) +
#'      ggtitle("Ankylosing Spondylitis Forest Plot",
#'              subtitle="Control Group Response Rate")
#'
#' # the defaults are set with
#' color_scheme_set("blue")
#' theme_set(theme_default(base_size=12))
#'
#' @export
forest_plot <- function(x,
                        prob=0.95,
                        est = c("both", "MAP", "Mean", "none"),
                        model = c("stratified", "both", "meta"),
                        point_est = c("median", "mean"),
                        size=1.25,
                        alpha=0.5) {

    assert_number(prob, lower=0, upper=1)
    assert_that(inherits(x, "gMAP"))
    assert_that(x$has_intercept)
    est <- match.arg(est)
    low <- (1-prob)/2
    up <- 1-low
    strat <- as.data.frame(x$est_strat(1-prob))
    strat <- cbind(strat[1:2], median=strat$mean, strat[3:4])
    names(strat)[3:4] <- c("low", "up")
    fit <- as.data.frame(fitted(x, type="response", probs=c(0.5, low, up)))

    est <- match.arg(est)
    model <- match.arg(model)
    point_est <- match.arg(point_est)

    if(est   == "both") est   <- c("MAP", "Mean")
    if(model == "both") model <- c("stratified", "meta")

    pred_est <- as.data.frame(do.call(rbind, summary(x, probs=c(0.5, low, up), type="response")[c("theta.pred", "theta")]))
    pred_est <- transform(pred_est,  study=c("MAP", "Mean") , model="meta")
    pred_est <- pred_est[c("MAP", "Mean") %in% est,]

    names(pred_est)[1:5] <- names(strat) <- names(fit) <- c("mean", "sem", "median", "low", "up")
    comb <- rbind(if("stratified" %in% model) transform(strat, study=rownames(strat), model="stratified"),
                  if("meta"       %in% model) transform(fit,   study=rownames(strat), model="meta"),
                  pred_est
                  )
    comb <- within(comb, { model <- factor(model, levels=c("meta", "stratified"))
                           study <- factor(study, levels=rev(c(rownames(strat), "Mean", "MAP"))) })

    opts <- list(position=position_dodge(width=0.3), size=size)

    xlab_str <- switch(x$family$family,
                       gaussian="Response",
                       binomial="Response Rate",
                       poisson="Counting Rate")

    graph <- ggplot(comb, aes_string(x="study", y=point_est, ymin="low", ymax="up", linetype="model", color="model"))

    if(any(c("MAP", "Mean") %in% est)) {
        ref_line <- est[est %in% c("Mean", "MAP")][1]
        ref_data <- subset(pred_est, study == ref_line)
        no_ref <- sum(est %in% c("Mean", "MAP"))
        graph <- graph + geom_rect(ymin=-Inf, ymax=Inf, xmin=0, xmax=no_ref + 0.5,
                                   fill=get_color("l"),
                                   color=get_color("l"), show.legend=FALSE) +
            geom_hline(yintercept=ref_data[1,point_est],
                       color=get_color("mh"),
                       alpha=alpha,
                       size=size)
    }

    graph <- graph +
        scale_color_manual("Model", values=get_color(c("mh", "m"))) +
        do.call(geom_pointrange, opts) +
        ylab(xlab_str) +
        scale_linetype_discrete("Model") +
        theme(axis.line.y=element_blank(), axis.ticks.y=element_blank()) +
        bayesplot::bayesplot_theme_get() +
        bayesplot::xaxis_title(FALSE) +
        coord_flip() +
        bayesplot::legend_none()

    graph
}
