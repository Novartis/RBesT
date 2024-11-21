#' internal function used for integration of densities which appears
#' to be much more stable from -Inf to +Inf in the logit space while
#' the density to be integrated recieves inputs from 0 to 1 such that
#' the inverse distribution function must be used. The integral solved
#' is int_x dmix(mix,x) integrand(x) where integrand must be given as
#' log and we integrate over the support of mix.
#'
#' integrate density in logit space and split by component such
#' that the quantile function of each component is used. This
#' ensures that the R implementation of the quantile function is
#' always used.
#'
#' @param log_integrand function to integrate over which must return the log(f)
#' @param mix density over which to integrate
#' @param Lplower logit of lower cumulative density
#' @param Lpupper logit of upper cumulative density
#'
#' @keywords internal
integrate_density_log <- function(log_integrand, mix, Lplower=-Inf, Lpupper=Inf, eps=getOption("RBesT.integrate_prob_eps", 1E-6)) {
    .integrand_comp_logit <- function(mix_comp) {
        function(l) {
            u   <- inv_logit(l)
            lp  <- log_inv_logit(l)
            lnp <- log_inv_logit(-l)
            exp(lp + lnp + log_integrand(qmix(mix_comp, u)))
        }
    }

    Nc <- ncol(mix)

    ## integrate by component of mix separatley to increase precision
    ## when the density is not 0 at the boundaries integration, then
    ## the integration is performed on the natural scale. The check
    ## for that is done on the identity scale to avoid numerical
    ## issues.
    lower  <- inv_logit(Lplower)
    upper  <- inv_logit(Lpupper)
    return(sum(vapply(1:Nc, function(comp) {
        mix_comp  <- mix[[comp, rescale=TRUE]]
        fn_integrand_comp_logit <- .integrand_comp_logit(mix_comp)
        if (all(!is.na(fn_integrand_comp_logit(c(Lplower, Lpupper))))) {
            return(.integrate(fn_integrand_comp_logit, Lplower, Lpupper))
        }
        lower_comp  <- ifelse(Lplower==-Inf, qmix(mix_comp, eps), qmix(mix_comp, lower))
        upper_comp  <- ifelse(Lpupper==Inf, qmix(mix_comp, 1-eps), qmix(mix_comp, upper))
        return(.integrate(function(x) exp(log_integrand(x) + dmix(mix_comp, x, log=TRUE)), lower_comp, upper_comp))
    }, c(0.1)) * mix[1,]))
}

integrate_density <- function(integrand, mix, Lplower=-Inf, Lpupper=Inf, eps=getOption("RBesT.integrate_prob_eps", 1E-6)) {
    .integrand_comp_logit <- function(mix_comp) {
        function(l) {
            u   <- inv_logit(l)
            lp  <- log_inv_logit(l)
            lnp <- log_inv_logit(-l)
            exp(lp + lnp) * integrand(qmix(mix_comp, u))
        }
    }
    Nc <- ncol(mix)

    lower  <- inv_logit(Lplower)
    upper  <- inv_logit(Lpupper)
    
    return(sum(vapply(1:Nc, function(comp) {
        mix_comp  <- mix[[comp, rescale=TRUE]]
        ## ensure that the integrand is defined at the boundaries...
        fn_integrand_comp_logit <- .integrand_comp_logit(mix_comp)
        if (all(!is.na(fn_integrand_comp_logit(c(Lplower, Lpupper))))) {
            return(.integrate(fn_integrand_comp_logit, Lplower, Lpupper))
        }
        ## ... otherwise we avoid the boundaries by eps prob density:
        lower_comp  <- ifelse(Lplower==-Inf, qmix(mix_comp, eps), qmix(mix_comp, lower))
        upper_comp  <- ifelse(Lpupper==Inf, qmix(mix_comp, 1-eps), qmix(mix_comp, upper))
        return(.integrate(function(x) integrand(x) * dmix(mix_comp, x), lower_comp, upper_comp))
    }, c(0.1)) * mix[1,]))
}

.integrate <- function(integrand, lower, upper) {
    integrate_args_user <- getOption("RBesT.integrate_args", list())
    args <- modifyList(list(lower=lower, upper=upper,
                            rel.tol=.Machine$double.eps^0.25,
                            abs.tol=.Machine$double.eps^0.25,
                            subdivisions=1000,
                            stop.on.error=TRUE),
                       integrate_args_user)

    integrate(integrand,
              lower=args$lower,
              upper=args$upper,
              rel.tol=args$rel.tol,
              abs.tol=args$abs.tol,
              subdivisions=args$subdivisions,
              stop.on.error=args$stop.on.error)$value
}
