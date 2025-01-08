#' ---
#' title: "Using RBesT to reproduce Schmidli et al. \"Robust MAP Priors\""
#' author: "Sebastian Weber"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Using RBesT to reproduce Schmidli et al. "Robust MAP Priors"}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteDepends{foreach}
#'   %\VignetteDepends{reshape2}
#' ---
#'
#' The following script demonstrates the ` RBesT` library to reproduce the
#' main results from Schmidli et al., Biometrics 70, 1024, 2014.
#'
#' The two main ideas of the paper are
#'
#' 1. use mixture priors to approximate accuratley numerical MCMC MAP
#'    priors
#'
#' 2. robustify informative MAP priors by adding a suitable
#'    non-informative component to the informative MAP prior
#'
#' As example an adaptive design for a binomial endpoint is considered:
#'
#' - Stage 1: mI in test treatment and nI in control (e.g., mI = 20, nI = 15);
#' - Stage 2: (m - mI ) in test treatment and max(n - ESSI, nmin) in control (e.g., nmin = 5).
#'
#' To render a report with all the figures, please run:
#'
#' - ` rmarkdown::render("robustMAP.R")`
#'
#' - ` rmarkdown::render("robustMAP.R", "pdf_document")` (for a pdf version)
#'

library(foreach)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RBesT)
library(knitr)
theme_set(theme_bw())
knitr::opts_chunk$set(
  fig.width = 7,
  fig.height = 4
)

#'
#' ## Operating Characteristics, Table 1
#'

#' the different priors
map <- list()
map$beta <- mixbeta(c(1.0, 4, 16))
map$mix90 <- mixbeta(c(0.9, 4, 16), c(0.1, 1, 1))
## map$mix70 <- mixbeta(c(0.7, 4, 16), c(0.3, 1, 1))
map$mix50 <- mixbeta(c(0.5, 4, 16), c(0.5, 1, 1))
## map$mix30 <- mixbeta(c(0.3, 4, 16), c(0.7, 1, 1))
## map$mix10 <- mixbeta(c(0.1, 4, 16), c(0.9, 1, 1))
map$unif <- mixbeta(c(1.0, 1, 1))

unif <- map$unif

## for the adaptive design the calculation is a bit involved as we
## have to calculate all possible ctl sample sizes which is determined
## by the ESS at the intermediate

OC_adaptBinary2S <- function(N1, Ntarget, Nmin, M, ctl.prior, treat.prior, pc, pt, decision) {
  ## calculate the different possible ESS values which we get after
  ## stage1
  r1 <- 0:N1
  ESSstage1 <- vector("double", N1 + 1)
  for (r in r1) {
    ESSstage1[r + 1] <- round(ess(postmix(ctl.prior, r = r, n = N1), method = "morita"))
  }

  ## number of patients enrolled in stage 2
  N2 <- pmax(Ntarget - ESSstage1, Nmin)

  ## total number of patients enrolled
  N <- N1 + N2

  P <- try(data.frame(pc = pc, pt = pt))
  if (inherits(P, "try-error")) {
    stop("pc and pt need to be of same size")
  }

  ## calculate for each scenario and sample size of the control the
  ## power
  power_all <- matrix(0, N1 + 1, nrow(P))
  for (r in r1) {
    ## power_all[r+1,] <- OC_binary2S(N[r+1], M, ctl.prior, treat.prior, P$pc, P$pt, crit)
    design_calc <- oc2S(treat.prior, ctl.prior, M, N[r + 1], decision)
    power_all[r + 1, ] <- design_calc(P$pt, P$pc)
  }

  ## finally take the mean with the respective weight which corresponds
  ## to the weight how the respective sample size occur
  w <- sapply(P$pc, function(p) dbinom(r1, N1, p))

  data.frame(power = colSums(power_all * w), samp = colSums(w * N))
}


## minimum of patients enrolled in stage 2
Nmin <- 5
## number of patients enrolled to control in stage 1
N1 <- 15
## target number of patients overall
Ntarget <- 40
## target number of patients in treatment group
M <- 40

## decision function: P(x1 - x2 > 0) > 0.975
dec <- decision2S(0.975, 0, lower.tail = FALSE)

cases <- expand.grid(prior = names(map), pc = seq(0.1, 0.6, by = 0.1), delta = c(0, 0.3))

## the mixture cases have a varying ess at the interim and need the
## adaptive function ...
cases.mix <- grep("mix", names(map), value = TRUE)
cases.fix <- grep("mix", names(map), value = TRUE, invert = TRUE)

resMix <- foreach(i = cases.mix, .combine = rbind) %do% {
  design <- subset(cases, prior == i)
  cbind(design, OC_adaptBinary2S(N1, Ntarget, Nmin, M, map[[i]], unif, design$pc, design$pc + design$delta, dec))
}

## ... for the non-mixture priors the ESS is fixed at the intermediate
## step such that the much faster oc2S can be used directly
resFix <- foreach(i = cases.fix, .combine = rbind) %do% {
  design <- subset(cases, prior == i)
  prior <- map[[i]]
  Nc <- Ntarget - round(ess(prior, method = "morita"))
  design_calc <- oc2S(unif, prior, M, Nc, dec)
  cbind(design, power = design_calc(design$pc + design$delta, design$pc), samp = Nc)
}

powerTable <- rbind(resMix, resFix)

P <- expand.grid(pc = c(0.2, 0.3, 0.4, 0.5), pt = seq(0.05, 0.95, by = 0.025))

powerMix <- foreach(i = cases.mix, .combine = rbind) %do% {
  cbind(P, prior = i, OC_adaptBinary2S(N1, Ntarget, Nmin, M, map[[i]], unif, P$pc, P$pt, dec))
}

powerFix <- foreach(i = cases.fix, .combine = rbind) %do% {
  prior <- map[[i]]
  Nc <- Ntarget - round(ess(prior, method = "morita"))
  design_calc <- oc2S(unif, prior, M, Nc, dec)
  cbind(P, prior = i, power = design_calc(P$pt, P$pc), samp = Nc)
  ## cbind(P, prior=i, power=OC_binary2S(Nc, M, prior, unif, P$pc, P$pt, dec), samp=Nc)
}

power <- rbind(powerMix, powerFix)


ocAdapt <- powerTable[, -ncol(powerTable)] %>%
  unite(case, delta, prior) %>%
  transform(power = 100 * power) %>%
  spread(case, power)
ocAdaptSamp <- powerTable[, -(ncol(powerTable) - 1)] %>%
  unite(case, delta, prior) %>%
  spread(case, samp)

kable(ocAdapt, digits = 1, caption = "Type I error and power")

kable(ocAdaptSamp, digits = 1, caption = "Sample size")

#'
#' ## Additional power Figure under varying pc
#'

ggplot(power, aes(pt - pc, power, colour = prior)) +
  geom_line() +
  facet_wrap(~pc) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(-0.8, 0.8, by = 0.2)) +
  coord_cartesian(xlim = c(-0.15, 0.5)) +
  geom_hline(yintercept = 0.025, linetype = 2) +
  geom_hline(yintercept = 0.8, linetype = 2) +
  ggtitle("Prob. for alternative for different pc")

#'
#' ## Bias and rMSE, Figure 1
#'
#' Reproduction of Fig. 1 in Robust MAP Prior paper.
#'


plot(map$beta, prob = 1)
plot(map$mix50)

#' The bias and rMSE calculations are slightly involved as the sample
#' size depends on the first stage.

est <- foreach(case = names(map), .combine = rbind) %do% {
  ## prior to consider
  prior <- map[[case]]

  ## calculate the different possible ESS values which we get after
  ## stage1
  r1 <- 0:N1
  ESSstage1 <- c()
  for (r in r1) {
    ESSstage1 <- c(ESSstage1, round(ess(postmix(prior, r = r, n = N1), method = "morita", loc = "mode")))
  }

  ## number of patients enrolled in stage 2
  N2 <- pmax(Ntarget - ESSstage1, 5)

  ## total number of patients enrolled
  N <- N1 + N2

  ## we need the maximal possible number of patients
  Nmax <- max(N)

  ## calculate for each i = 0 to N1 possible responders in stage one
  ## the posterior when observing 0 to N[i] in total. Calculate for
  ## each scenario outcome E(p) and E(p^2)
  m <- matrix(0, N1 + 1, Nmax + 1)
  m2 <- matrix(0, N1 + 1, Nmax + 1)
  for (i in seq_along(N)) {
    n <- N[i]
    for (r in 0:n) {
      res <- summary(postmix(prior, r = r, n = n))[c("mean", "sd")]
      m[i, r + 1] <- res["mean"]
      m2[i, r + 1] <- res["sd"]^2 + m[i, r + 1]^2
    }
  }

  ## now collect the terms correctly weighted for each assumed true rate
  bias <- rMSE <- c()
  pt <- seq(0, 1, length = 101)
  for (p in pt) {
    ## weight for each possible N at stage 1
    wnp <- dbinom(0:N1, N1, p)

    ## E(p) and E(p^2) for each possible N at stage 1
    Mnm <- rep(0, N1 + 1)
    Mnm2 <- rep(0, N1 + 1)

    ## for a given weight at stage 1....
    for (i in seq(N1 + 1)) {
      n <- N[i]
      ## weights of possible outcomes when having n draws in
      ## stage1, we go up to Nmax+1 to get a vector of correct
      ## length; all entries above n are set to 0 from dbinom as
      ## expected as we can never observe more counts than the
      ## number of trials...
      wp <- dbinom(0:Nmax, n, p)

      Mnm[i] <- sum(m[i, ] * wp)
      Mnm2[i] <- sum(m2[i, ] * wp)
    }

    ## ... which we average over possible outcomes in stage 1
    Mm <- sum(wnp * Mnm)
    Mm2 <- sum(wnp * Mnm2)

    bias <- c(bias, (Mm - p))
    rMSE <- c(rMSE, sqrt(Mm2 - 2 * p * Mm + p^2))
  }
  data.frame(p = pt, bias = bias, rMSE = rMSE, prior = case)
}


ggplot(est, aes(p, 100 * bias, colour = prior)) +
  geom_line() +
  ggtitle("Bias")
ggplot(est, aes(p, 100 * rMSE, colour = prior)) +
  geom_line() +
  ggtitle("rMSE")


#'
#' ##  Ulcerative colitis example
#'
#' Clinical example to exemplify the methodology.
#'

## set seed to guarantee exact reproducible results
set.seed(25445)

map <- gMAP(cbind(r, n - r) ~ 1 | study,
  family = binomial,
  data = colitis,
  tau.dist = "HalfNormal",
  beta.prior = 2,
  tau.prior = 1
)

map_auto <- automixfit(map)

## advanced: look at AIC of all fitted models
sapply(attr(map_auto, "models"), AIC)

print(map_auto)

## use best fitting mixture model as prior
prior <- map_auto

pl <- plot(prior)
pl$mix + ggtitle("MAP prior for ulcerative colitis")


#'
#' Colitis MAPs from paper for further figures.
#'

mapCol <- list(
  one = mixbeta(c(1, 2.3, 16)),
  two = mixbeta(c(0.77, 6.2, 50.8), c(1 - 0.77, 1.0, 4.7)),
  three = mixbeta(c(0.53, 2.5, 19.1), c(0.38, 14.6, 120.2), c(0.08, 0.9, 2.9))
)
mapCol <- c(mapCol, list(
  twoRob = robustify(mapCol$two, weight = 0.1, mean = 1 / 2),
  threeRob = robustify(mapCol$three, weight = 0.1, mean = 1 / 2)
))

#'
#' Posterior for different remission rates, Figure 3
#'


N <- 20
post <- foreach(prior = names(mapCol), .combine = rbind) %do% {
  res <- data.frame(mean = rep(NA, N + 1), sd = 0, r = 0:N)
  for (r in 0:N) {
    res[r + 1, 1:2] <- summary(postmix(mapCol[[prior]], r = r, n = N))[c("mean", "sd")]
  }
  res$prior <- prior
  res
}

ggplot(post, aes(r, mean, colour = prior, shape = prior)) +
  geom_point() +
  geom_abline(slope = 1 / 20)
ggplot(post, aes(r, sd, colour = prior, shape = prior)) +
  geom_point() +
  coord_cartesian(ylim = c(0, 0.17))

sessionInfo()
