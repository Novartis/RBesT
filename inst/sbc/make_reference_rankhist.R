#! /usr/bin/env Rscript

start_time  <- Sys.time()

here::i_am("inst/sbc/make_reference_rankhist.R")
library(here)

setwd(here())
system("make binary")
setwd(here("inst", "sbc"))

pkg <- c("assertthat", "rstan", "mvtnorm", "checkmate", "Formula", "abind", "dplyr", "tidyr", "here", "bayesplot")
sapply(pkg, require, character.only=TRUE)


library(clustermq)
library(data.table)
library(knitr)
sbc_tools <- new.env()
source(here("inst", "sbc", "sbc_tools.R"), local=sbc_tools)
set.seed(453453)

scheduler <- getOption("clustermq.scheduler")

if(is.null(scheduler)) {
    ## in this case we enable the multiprocess option to leverage local CPUs
    options(clustermq.scheduler="multiprocess")
}

scheduler <- getOption("clustermq.scheduler")

##options(clustermq.scheduler="LOCAL")
n_jobs <- 1
if(scheduler == "multiprocess") {
    ## on a local machine we only use as many CPUs as available
    n_jobs  <- as.numeric(system2("nproc", stdout=TRUE))
}
if(scheduler %in% c("LSF", "SGE", "SLURM", "PBS", "Torque")) {
    ## on a queinging enabled backend, we use a lot more parallelism
    n_jobs  <- 200
}

cat("Using clustermq backend", scheduler, "with", n_jobs, "concurrent jobs.\n")

#' Evaluate dense and sparse data-scenario

#' - Dense: 10 trials with 40 entries each
#' - Sparse: 3 trials with 40 entries each
base_scenarios <- list(dense=list(group=rep(1:10, each=40)),
                       sparse=list(group=rep(1:3, each=40)))

## family, mean_mu, sd_mu, sd_tau, samp_sd
cases <- data.frame(
    family=c("binomial", "gaussian", "poisson"),
    mean_mu=c(-1, 0, 0),
    sd_mu=c(1),
    sd_tau=c(rep(0.5, 3), rep(1, 3)),
    samp_sd=c(1),
    stringsAsFactors=FALSE)

## replications to use
S <- 1E4

scenarios <- merge(expand.grid(repl=1:S, data_scenario=c("dense", "sparse"), stringsAsFactors=FALSE), cases, by=NULL)
##scenarios <- merge(expand.grid(repl=1:S, data_scenario=c("sparse"), stringsAsFactors=FALSE), cases, by=NULL)

scenarios <- cbind(job.id=1:nrow(scenarios), scenarios)

num_simulations <- nrow(scenarios)

cat("Total number of jobs to dispatch:", num_simulations, "\n")

RNGkind("L'Ecuyer-CMRG")
set.seed(56969)
rng_seeds <- sbc_tools$setup_lecuyer_seeds(.Random.seed, num_simulations)

sim_result <- Q_rows(scenarios, sbc_tools$run_sbc_case, const=list(base_scenarios=base_scenarios, seeds=rng_seeds), export=as.list(sbc_tools), n_jobs=n_jobs, pkgs=pkg)

assert_that(num_simulations == length(sim_result), msg="Check if all simulations were processed.")

calibration_data <- merge(scenarios, bind_rows(sim_result), by="job.id")

## convert to data.table
setDT(calibration_data)

## collect sampler diagnostics
sampler_diagnostics <- calibration_data %>%
    group_by(family, data_scenario, sd_tau) %>%
    summarize(N=n(),
              total_divergent=sum(n_divergent),
              total_divergent_sim_fraction=mean(n_divergent>0),
              min_ess=min(min_Neff),
              max_Rhat=max(max_Rhat),
              total_large_Rhat=sum(max_Rhat > 1.2),
              min_lp_ess_bulk=min(lp_ess_bulk),
              min_lp_ess_tail=min(lp_ess_tail))


cat("\nSampler diagnostics:\n\n")
kable(sampler_diagnostics, digits=3)
cat("\n")

if(sum(sampler_diagnostics$total_divergent) != 0) {
    warning("There were some divergent transitions!")
}
if(any(sampler_diagnostics$max_Rhat > 1.2) ) {
    warning("There were some parameters with large Rhat!")
}

#' Bin raw data as used in the analysis.
scale64  <- sbc_tools$scale_ranks(1024, 2^4)
B <- 1024L / 2^4
calibration_data_binned <- calibration_data[, scale64(.SD), by=c("data_scenario", "family", "sd_tau")]

#' Save as data.frame to avoid data.table dependency.
calibration_data <- as.data.frame(calibration_data)
calibration_data_binned <- as.data.frame(calibration_data_binned)

#' Further identification and verification data of run
git_hash <- system2("git", c("rev-parse", "HEAD"), stdout=TRUE)
created <- Sys.time()
created_str <- format(created, "%F %T %Z", tz="UTC")

calibration <- list(## raw=calibration_data, ## stop storing raw results, which are not needed for SBC reports
                    data=calibration_data_binned,
                    sampler_diagnostics = sampler_diagnostics,
                    S=S,
                    B=B,
                    git_hash=git_hash,
                    created=created)

saveRDS(calibration, file="calibration.rds")
saveRDS(calibration_data, file="calibration_data.rds")

library(tools)
md5 <- md5sum("calibration.rds")
cat(paste0("Created:  ", created_str, "\ngit hash: ", git_hash, "\nMD5:      ", md5, "\n"), file="calibration.md5")

#'
#' Summarize execution time
#'
job_report <- calibration_data[c("job.id", "time.running", names(scenarios))]
setDT(job_report)
job_report$time.running <-  job_report$time.running / 60 ## convert to minutes

runtime_by_problem_family  <- job_report %>%
    group_by(family, data_scenario) %>%
    summarize(total=sum(time.running), mean=mean(time.running), max=max(time.running))

runtime_by_problem  <- job_report %>%
    group_by(data_scenario) %>%
    summarize(total=sum(time.running), mean=mean(time.running), max=max(time.running))

runtime  <- job_report %>%
    group_by(family) %>%
    summarize(total=sum(time.running), mean=mean(time.running), max=max(time.running))

cat("Summary on job runtime on cluster:\n\n")

cat("\nRuntime by problem and family:\n")
kable(runtime_by_problem_family, digits=2)

cat("\nRuntime by family:\n")
kable(runtime, digits=2)

cat("\nRuntime by problem:\n")
kable(runtime_by_problem, digits=2)

end_time <- Sys.time()

total_runtime <- difftime(end_time, start_time)
units(total_runtime) <- "mins"

cat("\n\nTotal runtime (min):", as.numeric(total_runtime), "\n\n\n")

#' Session info
sessionInfo()
