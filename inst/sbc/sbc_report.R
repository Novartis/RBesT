#' ---
#' title: Simulation based calibration for RBesT
#' author: "Sebastian Weber"
#' date: "`r date()`"
#' output: html_vignette
#' params:
#'   include_plots: FALSE
#' vignette: >
#'    %\VignetteIndexEntry{Simulation based calibration for RBesT}
#'    %\VignetteEngine{knitr::rmarkdown}
#'    %\VignetteEncoding{UTF-8}
#' ---
#'
#+ include=FALSE
library(knitr)
library(tools)
library(assertthat)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
theme_set(theme_bw())
here::i_am("inst/sbc/sbc_report.R")
library(here)
source(here("inst", "sbc", "sbc_tools.R"))
library(rstan)
library(purrr)

knitr::opts_chunk$set(
  fig.width = 1.62 * 4,
  fig.height = 4,
  cache = FALSE,
  echo = FALSE
)
#'
#' This report documents the results of a simulation based calibration
#' (SBC) run for `RBesT`. The calibration data will be generated
#' whenever relevant changes to the `gMAP` function were made. The
#' calibration runs are performed for typical use cases of
#' `gMAP`. These include the three likelihoods (binomial, gaussian &
#' Poisson), a sparse ($2$ trials) and dense ($10$ trials) data
#' situation and finally a run with a very/less conservative prior
#' choice for between-trial heterogeniety parameter.
#'
#' The calibration data presented here has been generated at and with
#' the `RBesT` git version as:
cat(readLines(here("inst", "sbc", "calibration.md5")), sep = "\n")
#'
#' The MD5 hash of the calibration data file presented here must match
#' the above listed MD5:
md5sum(here("inst", "sbc", "calibration.rds"))
#'
#' # Introduction
#'
#' Simulation based calibration (SBC) is a necessary condition which
#' must be met for any Bayesian analysis with proper priors. The
#' details are presented in Talts, et. al (see
#' https://arxiv.org/abs/1804.06788).
#'
#' Self-consistency of any Bayesian analysis with a proper prior:
#'
#' $$ p(\theta) = \iint \mbox{d}\tilde{y} \, \mbox{d}\tilde{\theta} \, p(\theta|\tilde{y}) \, p(\tilde{y}|\tilde{\theta}) \, p(\tilde{\theta}) $$
#' $$ \Leftrightarrow p(\theta) = \iint \mbox{d}\tilde{y} \, \mbox{d}\tilde{\theta} \, p(\theta,\tilde{y},\tilde{\theta}) $$
#'
#' SBC procedure:
#'
#' Repeat $s=1, ..., S$ times:
#'
#' 1. Sample from the prior $$\tilde{\theta} \sim p(\theta)$$
#'
#' 2. Sample fake data $$\tilde{y} \sim p(y|\tilde{\theta})$$
#'
#' 3. Obtain $L$ posterior samples $$\{\theta_1, ..., \theta_L\} \sim p(\tilde{\theta}|\tilde{y})$$
#'
#' 4. Calculate the *rank* $r_s$ of the prior draw $\tilde{\theta}$ wrt to
#' the posterior sample $\{\theta_1, ..., \theta_L\} \sim p(\tilde{\theta}|\tilde{y})$ which falls into the range $[0,L]$
#' out of the possible $L+1$ ranks. The rank is calculated as
#' $$r_s = \sum_{l=1}^L \mathbb{I}[ \theta_l < \tilde{\theta}]$$
#'
#' The $S$ ranks then form a uniform $0-1$ density and the count in
#' each bin has a binomial distribution with probability of
#' $$p(r \in \mbox{Any Bin}) =\frac{(L+1)}{S}.$$
#'
#' ## Hierarchical intercept only (random-effects intercept) model for binomial, gaussian and Poisson likelihood
#'
#' Likelihood:
#'
#' - Binary $$y_i|\theta_{j} \sim \mbox{Bernoulli}(\theta_j), $$ $$g(\theta) = \mbox{logit}(\theta)$$
#' - Normal $$y_i|\theta_{j} \sim \mbox{Normal}(\theta_j, \sigma^2), $$ $$g(\theta) = \theta$$
#' - Poisson $$y_i|\theta_{j} \sim \mbox{Poisson}(\theta_j), $$ $$g(\theta) = \log(\theta)$$
#'
#' Hierarchical prior:
#'
#' $$ g(\theta_j)|\mu,\tau \sim \mbox{Normal}(\mu, \tau^2)$$
#'
#' $$\mu \sim \mbox{Normal}(m_\mu, s^2_\mu)$$
#' $$\tau \sim \mbox{Normal}^+(0, s^2_\tau)$$
#'
#' The fake data simulation function returns for binomial and Poisson
#' data the sum of the responses while for normal the mean summary is
#' used. Please refer to the `sbc_tools.R` and
#' `make_reference_rankhist.R` R programs for the implementation
#' details.
#'
#' The reference runs are created with $L=1023$ posterior draws for
#' each replication and a total of $S=10^4$ replications are run per
#' case. For the evaluation here the results are reduced to
#' $B=L'+1=64$ bins to ensure a sufficiently large sample size per
#' bin.
#'

calibration <- readRDS(here("inst", "sbc", "calibration.rds"))
include_plots <- TRUE
if ("params" %in% ls()) {
  include_plots <- params$include_plots
}

# The summary function we use here scales down the $L+1=1024$ bins to
# smaller number of rank bins. This improves the number of counts
# expected per rank bin ($S/(L+1)$) and is thus more robust in terms
# of large number laws. We choose $L=1023$ samples from the posterior
# such that we have $1024 = 2^10$ bins for the ranks. Thus any power
# of $2$ can be used to scale down the number of bins.

plot_binned <- function(count, rank, group) {
  S <- sum(count[group == group[1]])
  num_ranks <- length(rank[group == group[1]])
  c95 <- qbinom(c(0.025, 0.5, 0.975), S, 1 / num_ranks)
  dd <- arrange(data.frame(count = count, rank = rank, group = group, stringsAsFactors = FALSE), group, rank) %>%
    group_by(group) %>%
    mutate(ecdf = cumsum(count) / S, ecdf_ref = (rank + 1) / (num_ranks)) %>%
    separate(group, into = c("g1", "g2"), sep = "/")
  pl <- list()
  pl[["hist"]] <- ggplot(dd, aes(rank, count)) +
    facet_grid(g1 ~ g2) +
    geom_col() +
    geom_hline(yintercept = c95[c(1, 3)], linetype = I(2)) +
    geom_hline(yintercept = c95[c(2)], linetype = I(3))
  pl[["ecdf_diff"]] <- ggplot(dd, aes(rank, ecdf - ecdf_ref)) +
    facet_grid(g1 ~ g2) +
    geom_step() +
    geom_hline(yintercept = 0, linetype = I(3))
  pl
}


B <- calibration$B
S <- calibration$S

calibration_binned <- calibration$data %>%
  unite(family, sd_tau, col = "group", sep = "/") %>%
  group_by(data_scenario, group)

calibration_binned <- calibration_binned %>%
  gather(starts_with("count"), key = "parameter", value = "count") %>%
  mutate(parameter = sub("^count.", "", parameter))

## filter out cases where count == 0 for all entries of the parameters
## (happens due the way data is processed for the sparse cases)

calibration_binned <- calibration_binned %>%
  group_by(data_scenario, group, parameter) %>%
  mutate(all_zero = all(count == 0)) %>%
  ungroup() %>%
  subset(!all_zero) %>%
  mutate(all_zero = NULL)

calibration_dense <- subset(calibration_binned, data_scenario == "dense")
calibration_sparse <- subset(calibration_binned, data_scenario == "sparse")

pl_dense <- calibration_dense %>%
  split(.$parameter) %>%
  map(~ plot_binned(.$count, .$rank, .$group))

pl_sparse <- calibration_sparse %>%
  split(.$parameter) %>%
  map(~ plot_binned(.$count, .$rank, .$group))

#' # SBC results
#'
#' ## Sampler Diagnostics Overview
#'

kable(calibration$sampler_diagnostics, digits = 3)

#'
#' Note: Large Rhat is defined as exceeding 1.2.
#'

#'
#' # Summary Statistics
#'
#' ## $\chi^2$ Statistic, $\mu$
#'

chisq <- calibration_binned %>%
  group_by(data_scenario, group, parameter) %>%
  group_map(~ cbind(case = .y, tidy(chisq.test(.$count))[, c(1, 3, 2)])) %>%
  bind_rows() %>%
  rename(df = parameter, data_scenario = case.data_scenario, group = case.group, parameter = case.parameter) %>%
  separate(group, into = c("likelihood", "sd_tau"), sep = "/")

kable(subset(chisq, parameter == "mu"), digits = 3)

#'
#' ## $\chi^2$ Statistic, $\tau$
#'

kable(subset(chisq, parameter == "tau"), digits = 3)

#'
#' ## $\chi^2$ Statistic, group estimates $\theta$
#'

kable(subset(chisq, parameter != "tau" & parameter != "mu"), digits = 3)

#+ results="asis", include=include_plots, eval=include_plots
spin_child("sbc_report_plots.R")

#'
#' ## Session Info
#'
sessionInfo()
