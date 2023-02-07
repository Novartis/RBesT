#' Diagnostic plots for gMAP analyses
#' 
#' @param x \code{\link{gMAP}} object
#' @param size Controls line sizes of traceplots and forest plot.
#' @param ... Ignored.
#'
#' @details Creates MCMC diagnostics and a forest plot (including
#' model estimates) for a \code{\link{gMAP}} analysis. For a
#' customized forest plot, please use the dedicated function
#' \code{\link{forest_plot}}.
#'
#' @template plot-help
#' 
#' @return The function returns a list of \code{\link{ggplot}}
#' objects.
#'
#' @method plot gMAP
#' @export
plot.gMAP <- function(x, size=NULL, ...) {
    pl <- list()

    draws_all <- rstan::extract(x$fit, permuted=FALSE, inc_warmup=TRUE)
    thin <- attr(x$fit, "sim")$thin
    n_warmup <- floor(attr(x$fit, "sim")$warmup / thin)
    draws <- draws_all[- (1:n_warmup),,]
    nuts_diag <- bayesplot::nuts_params(x$fit, inc_warmup=FALSE)

    ## by default we return a small set of plots only
    plot_verbose <- getOption("RBesT.verbose", FALSE)
    
    div_opts <- list()
    if(sum(nuts_diag$Value[nuts_diag$Parameter=="divergent__"]) > 0) {
        div_opts$np <- bayesplot::nuts_params(x$fit, inc_warmup=TRUE)
        plot_verbose <- TRUE ## if any divergent transition happens,
                             ## then we plot verbose in any case
    }

    tau_pars <- grep("tau\\[", dimnames(draws)$parameters, value=TRUE)
    beta_pars <- grep("beta\\[", dimnames(draws)$parameters, value=TRUE)

    tau_log_trans <- as.list(rep("log", length(tau_pars)))
    names(tau_log_trans) <- tau_pars    

    if(plot_verbose) {
        ## traces are only shown if in verbose mode...
        pl$traceBeta <- do.call(bayesplot::mcmc_trace, c(list(x=draws_all, pars=beta_pars, size=size, n_warmup=n_warmup,
                                                   facet_args = list(labeller = ggplot2::label_parsed)), div_opts)) +
            ggtitle(expression(paste("Trace of Regression Coefficient ", beta))) +
                bayesplot::facet_text(length(beta_pars)!=1) + xlab("Iteration")
        pl$traceTau <- do.call(bayesplot::mcmc_trace, c(list(x=draws_all, pars=tau_pars, size=size, n_warmup=n_warmup,
                                                  facet_args = list(labeller = ggplot2::label_parsed)), div_opts)) +
            bayesplot::facet_text(length(tau_pars)!=1) +
                ggtitle(expression(paste("Trace of Heterogeneity Parameter ", tau))) +
                    xlab("Iteration")
        pl$traceLogTau <- do.call(bayesplot::mcmc_trace, c(list(x=draws_all, pars=tau_pars, size=size, n_warmup=n_warmup,
                                                     facet_args = list(labeller = ggplot2::label_parsed), transformations=tau_log_trans), div_opts)) +
            bayesplot::facet_text(length(tau_pars)!=1) +
                ggtitle(expression(paste("Trace of Heterogeneity Parameter ", tau, " on log-scale"))) +
                    xlab("Iteration")

        ## ... as well as auxilary model parameters
        pl$densityBeta <- bayesplot::mcmc_dens_overlay(x=draws, pars=beta_pars, facet_args = list(labeller = ggplot2::label_parsed, strip.position="bottom")) +
            ggtitle(expression(paste("Density of Regression Coefficient ", beta)))
        pl$densityTau <- bayesplot::mcmc_dens_overlay(x=draws, pars=tau_pars, facet_args = list(labeller = ggplot2::label_parsed, strip.position="bottom")) +
            ggtitle(expression(paste("Density of Heterogeneity Parameter ", tau)))
        pl$densityLogTau <- bayesplot::mcmc_dens_overlay(x=draws, pars=tau_pars, facet_args = list(labeller = ggplot2::label_parsed, strip.position="bottom"), transformations=tau_log_trans) +
            ggtitle(expression(paste("Density of Heterogeneity Parameter ", tau, " on log-scale")))
    }

    if(x$has_intercept) {
        pl$densityThetaStar <- bayesplot::mcmc_dens_overlay(x=draws, pars="theta_resp_pred") + xlab(expression(theta[symbol("*")])) + bayesplot::facet_text(FALSE) + ggtitle(expression(paste("Density of MAP Prior ", theta[symbol("*")])))
        pl$densityThetaStarLink <- bayesplot::mcmc_dens_overlay(x=draws, pars="theta_pred") +  xlab(expression(theta[symbol("*")])) + bayesplot::facet_text(FALSE) + ggtitle(expression(paste("Density of MAP Prior ", theta[symbol("*")], " (link scale)")))

        pl$forest_model <- forest_plot(x, model="both", size=if(is.null(size)) 1.25 else size)
    } else {
        message("No intercept defined.")
    }
    
    return(pl)
}
