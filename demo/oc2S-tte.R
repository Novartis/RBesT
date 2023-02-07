## 2-sample operating characeristics for time-to-event analysis
##
## Assumptions:
## - fixed sample size (aka events) per arm
## - simulation of accrual process
## - trial is censored once a critical number of events is reached
## - Poisson likelihood used => time-units must be sufficiently precise
## - censoring only due to analysis cut-off (administrative)
## - todo: add drop-out censoring
## - constant hazard over time
##

library(RBesT)
library(assertthat)
library(dplyr)
library(purrr)
library(parallel)
library(checkmate)

## recommended
options(mc.cores=detectCores(logical=FALSE))

######## Utility functions ########

## simulate the follow-up time for a given sample size when censoring
## with event_crit events. Positive time-points correspond to events,
## negative ones are censored cases.
sim_trial <- function(n_1, n_2, event_rate_1, event_rate_2, events_crit, accrual_rate) {
    n <- n_1 + n_2
    assert_that(events_crit <= n)
    accrual <- cumsum(c(0, rexp(n-1, accrual_rate)))
    idx_1 <- sample.int(n, n_1)
    idx_2  <- setdiff(1:n, idx_1)
    event_1 <- rexp(n_1, event_rate_1)
    event_2 <- rexp(n_2, event_rate_2)
    event <- rep(0, n)
    group  <- rep(1, n)
    group[idx_2] <- 2
    event[idx_1] <- event_1
    event[idx_2] <- event_2
    event_calendar  <- accrual + event
    o <- order(event_calendar)
    cens_calendar <- event_calendar[o[events_crit]]
    is_event <- event_calendar <= cens_calendar
    is_enrolled <- cens_calendar >= accrual
    follow_up <- ifelse(is_event, event_calendar - accrual, -1 * (cens_calendar - accrual) )
    follow_up[!is_enrolled] <- 0
    cbind(follow_up=follow_up, group=group)
}

## given the follow_up from one arm calculates the posterior from the
## given prior
analyze_arm  <- function(prior, follow_up) {
    assert_set_equal(likelihood(prior), "poisson")
    ## total time-units
    exposure_time <- sum(ceiling(abs(follow_up)))
    ## total events
    num_events <- sum(follow_up > 0)
    postmix(prior, m=num_events/exposure_time, n=exposure_time)
}

## analyzes a simulated trial. Returns 0 for failure or 1 for success
## depending on the decision function.
analyze_trial <- function(trial, prior_1, prior_2, decision) {
    is_grp_1  <- trial[,"group"] == 1
    post_1 <- analyze_arm(prior_1, trial[is_grp_1, "follow_up"])
    post_2 <- analyze_arm(prior_2, trial[!is_grp_1, "follow_up"])
    decision(post_1, post_2)
}

#' @param prior_1,prior_2 prior for each arm
#' @param n_1,n_2 sample size per arm
#' @param events_crit critical number of events when analysis cut-off censoring occurs
#' @param accrual_rate accrual rate of subjects
#' @param decision 2-sample decision function
#' @param num_sim number of replicates for each scenario
oc2S_tte <- function(prior_1, prior_2,
                     n_1, n_2,
                     events_crit,
                     accrual_rate,
                     decision,
                     num_sim=1E3) {

    fn <- function(event_rate_1, event_rate_2) {
        sim <- replicate(num_sim, sim_trial(n_1, n_2, event_rate_1, event_rate_2, events_crit, accrual_rate))
        mean(apply(sim, 3, partial(analyze_trial,
                                   prior_1=prior_1, prior_2=prior_2,
                                   decision=decision)))
    }
    function(theta1, theta2) {
        T <- try(data.frame(theta1 = theta1, theta2 = theta2, row.names=NULL))
        if (inherits(T, "try-error")) {
            stop("theta1 and theta2 need to be of same size")
        }
        mcmapply(fn, T$theta1, T$theta2)
    }
}


## Example with one historical study

## Historical data (assume constant hazard) given in month time-units
hist_data <- data.frame(events=14, exposure_m=220, study ="HISTDATA")

## changing time-units so that we have about 1 event per time-unit =>
## so we use years in this case

days_in_year <- 365.25
months_in_year <- 12
days_in_month <- days_in_year/months_in_year

hist_data <- hist_data %>%
    mutate(exposure_y=exposure_m/months_in_year)

set.seed(456747)
mc_map <- gMAP(events ~ 1 + offset(log(exposure_y)) | study,
               data = hist_data,
               tau.dist = "LogNormal",
               tau.prior = c(log(0.125), log(2)/1.96), ## assuming moderate heterogeniety
               beta.prior = 1,
               family = "poisson")


plot(mc_map)$forest_model

map <- mixfit(mc_map, Nc=2)

plot(map)$mix

## to get the ESS, we need to change the likelihood to exp
likelihood(map) <- "exp"
ess(map)

## change it back to a poisson likelihood as used during the OCs
likelihood(map) <- "poisson"


## OC setup

## decision criteria (trt - ctl)
success  <- decision2S(0.975, 0, lower.tail=TRUE, link="log")

## weakly-informative prior for treatment (encoding that the events
## are on year scale when using days as basic unit)
prior_noninf  <- mixgamma(c(1, 1/days_in_year, 1), param="mn", likelihood="exp")
prior_noninf

## change to poisson likelihood
likelihood(prior_noninf)  <- "poisson"

## change map from years units to days
map_d  <- map
map_d["b",]  <- map_d["b",] * days_in_year

## event rate in units of per day; 0.07 is per month
lambda_ctl  <- 0.07 / days_in_month
lambda_trt  <- 0.04 / days_in_month

## test 100 vs 50 patients and trial is stopped at 100 events total.
## assuming a 15 pts / month accrual rate

## note: for debugging set use
## options(mc.cores=1)

## not using a prior
design_noninf <- oc2S_tte(prior_noninf, prior_noninf,
                          100, 50,  ## sample size per arm
                          100,      ## critical number of events when trial is censored
                          15/days_in_month,  ## accrual in # of pts per day
                          success
                          )

oc  <- data.frame(theta1=c(lambda_trt, lambda_ctl, lambda_trt),
                  theta2=c(lambda_ctl, lambda_ctl, lambda_trt)) %>%
    mutate(noninf=design_noninf(theta1, theta2))


## with map prior on ctl
summary(map_d)
lambda_ctl

design_map <- oc2S_tte(prior_noninf, map_d,
                       100, 50,  ## sample size per arm
                       100,      ## critical number of events when trial is censored
                       15/days_in_month,  ## accrual in # of pts per day
                       success
                       )

oc <- oc %>%
    mutate(map=design_map(theta1, theta2))

oc
