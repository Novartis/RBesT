
## test gMAP results using SBC and with matching rstanarm models

suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(dplyr))
suppressWarnings(suppressPackageStartupMessages(library(rstanarm)))

eps <- 2e-1

probs <- seq(0.1, 0.9, by=0.1)

fast_sampling <- list(RBesT.MC.warmup=250, RBesT.MC.iter=500, RBesT.MC.chains=2)
std_sampling  <- list(RBesT.MC.warmup=NULL, RBesT.MC.iter=NULL, RBesT.MC.chains=NULL, RBesT.MC.control=NULL, RBesT.MC.ncp=NULL)
bad_sampling  <- list(RBesT.MC.chains=1, RBesT.MC.control=list(adapt_delta=0.75, stepsize=1), RBesT.MC.ncp=0)

## standardize posterior by median and IQR
get_std_quants <- function(sim, med, disp) {
    if(missing(med))
        med <- median(sim)
    if(missing(disp))
        disp <- IQR(sim)/1.34
    sim_std <- scale(sim, center=med, scale=disp)
    ref_quants <- quantile(sim_std, prob=probs)
    list(ref=ref_quants, med=med, disp=disp)
}


## used for now for comparisons to rstanarm (maybe drop as we have SBC
## now?)
cmp_reference <- function(best_gmap, OB_ref) {
    best_sim <- rstan::extract(best_gmap$fit, pars=c("beta", "tau", "theta_resp_pred"))
    for(n in names(best_sim)) {
        OB_case <- OB_ref[[n]]
        case <- modifyList(OB_case, list(ref=NULL, sim=best_sim[[n]]))
        best_std <- do.call(get_std_quants, case)
        res <- abs(best_std$ref - OB_case$ref)
        ##cat("Testing", n, ":", signif(res, 3), "\n")
        ##cat(n, " = ", paste(round(res, 4), collapse=", "), "\n")
        expect_true(all(res < eps))
    }
}

make_rstanarm_ref <- function(stanreg) {
    pred_arm <- stanreg$family$linkinv(posterior_linpred(stanreg, newdata=data.frame(study="MAP", e=1)))
    arm_post <- cbind(as.matrix(stanreg, pars=c("(Intercept)", "Sigma[study:(Intercept),(Intercept)]")), pred_arm)
    arm_post[,2] <- sqrt(arm_post[,2])
    colnames(arm_post) <- c("beta", "tau", "theta_resp_pred")
    arm_ref <- lapply(1:3, function(v) get_std_quants(arm_post[,v]))
    names(arm_ref) <- colnames(arm_post)
    arm_ref
}


## these examples need very high adapt_delta to avoid divergent
## transitions
do.call(options, std_sampling)
options(RBesT.MC.control=list(adapt_delta=0.999, stepsize=0.01))

test_that("gMAP meets SBC requirements wrt to a Chi-Square statistic.", {
    require(dplyr)
    require(tidyr)
    sbc_chisq_test <-  RBesT:::calibration_data %>%
        gather(count.mu, count.tau, key="param", value="count") %>%
        group_by(data_scenario, family, sd_tau, param) %>%
        do(as.data.frame(chisq.test(.$count)[c("statistic", "p.value")]))
    num_tests  <- nrow(sbc_chisq_test)
    num_failed <- sum(sbc_chisq_test$p.value < 0.05)
    pv <- pbinom(num_failed, num_tests, 0.05)
    expect_true( pv > 0.025 & pv < 0.975 )
}
)

test_that("gMAP meets SBC requirements per bin.", {
    require(dplyr)
    require(tidyr)
    B <- RBesT:::calibration_meta$B
    S  <- RBesT:::calibration_meta$S
    alpha  <- 0.2
    ptrue  <- 1/B
    crit_low  <- qbinom(alpha/2, S, ptrue)
    crit_high  <- qbinom(1-alpha/2, S, ptrue)
    sbc_binom_test <-  RBesT:::calibration_data %>%
        gather(count.mu, count.tau, key="param", value="count") %>%
        group_by(data_scenario, family, sd_tau, param) %>%
        summarise(crit=sum(count < crit_low | count > crit_high)) %>%
        mutate(pvalue=pbinom(crit, B, alpha), extreme=pvalue<0.025|pvalue>0.975)
    num_tests  <- nrow(sbc_binom_test)
    num_failed <- sum(sbc_binom_test$extreme)
    pv <- pbinom(num_failed, num_tests, 0.05)
    expect_true( pv > 0.025 & pv < 0.975 )
}
)

test_that("SBC data was up to date at package creation.", {
    calibration_datum  <- RBesT:::calibration_meta$created
    package_datum <- RBesT:::pkg_create_date
    delta <- difftime(package_datum, calibration_datum, units="weeks")
    expect_true(delta < 52./2.)
}
)

## match against respective rstanarm model
set.seed(92575)
rate <- round(-log(0.05)/2, 1)
test_that("gMAP matches RStanArm binomial family", {
    skip("RStanArm has issues loading since 2024-01-02 on CI/CD systems.")
              skip_on_cran()
              suppressWarnings( best_run <- gMAP(cbind(r, n-r) ~ 1 | study,
                               data=AS,
                               family=binomial,
                               tau.dist="Exp",
                               tau.prior=c(rate),
                               beta.prior=cbind(0, 2)
                               ) )
              suppressWarnings( out <- capture.output(rstanarm_run <- make_rstanarm_ref(
                  stan_glmer(cbind(r, n-r) ~ 1 + (1|study),
                             data=AS,
                             family=binomial,
                             refresh=0,
                             iter=4000,
                             warmup=1000,
                             adapt_delta=0.999,
                             seed=4356,
                             chains=4,
                             prior=normal(0,2,autoscale=FALSE),
                             prior_intercept=normal(0,2,autoscale=FALSE),
                             prior_covariance=decov(1, 1, 1, 1/rate)
                             )))
                  )
              cmp_reference(best_gmap=best_run, OB_ref=rstanarm_run)
          })

## the remaining tests do not rely on good sampling, hence speed it up
do.call(options, fast_sampling)

## add test case with a single data
test_that("gMAP processes single trial case", {
              suppressMessages(suppressWarnings(map1 <- gMAP(cbind(r, n-r) ~ 1,
                           data=colitis[1,],
                           family=binomial,
                           tau.dist="HalfNormal",
                           ## prior are choosen super-tight to avoid sampling trouble
                           tau.prior=c(0.5),
                           beta.prior=cbind(0, 1)
                           )))
              expect_true(nrow(fitted(map1)) == 1)
          }
          )

test_that("gMAP processes not continuously labeled studies", {
              suppressWarnings( out <- capture.output(map1 <- gMAP(cbind(r, n-r) ~ 1 | study, data=AS[-1,],
                                                 family=binomial, tau.dist="HalfNormal", tau.prior=0.5,
                                                 iter=100, warmup=50, chains=1, thin=1))
                               )
              expect_true(nrow(fitted(map1)) == nrow(AS) - 1)
          })

## set bad sampling parameters to trigger divergences in the next test
do.call(options, bad_sampling)

set.seed(23434)
test_that("gMAP reports divergences", {
              suppressMessages(suppressWarnings(mcmc_div <- gMAP(cbind(r, n-r) ~ 1 | study, data=AS[1,,drop=FALSE], family=binomial,
                                                                 tau.dist="Uniform", tau.prior=cbind(0, 1000),
                                                                 beta.prior=cbind(0,1E5),
                                                                 iter=1000, warmup=0, chains=1, thin=1, init=10)))
              sp <- rstan::get_sampler_params(mcmc_div$fit, inc_warmup=FALSE)[[1]]
              expect_true(sum(sp[,"divergent__"]) > 0)
          })

## set sampling back to standards
do.call(options, std_sampling)

test_that("gMAP handles extreme response rates", {
              n <- 5
              data1 <- data.frame(n=c(n,n,n,n),r=c(5,5,5,5), study=1)
              suppressWarnings( map1 <- gMAP(cbind(r, n-r) ~ 1 | study, family=binomial,
                                             data=data1, tau.dist="HalfNormal",
                                             tau.prior=2.0, beta.prior=2,
                                             warmup=100, iter=200, chains=1, thin=1) )
              expect_true(nrow(fitted(map1)) == 4)
              data2 <- data.frame(n=c(n,n,n,n),r=c(0,0,0,0), study=1)
              suppressWarnings( map2 <- gMAP(cbind(r, n-r) ~ 1 | study, family=binomial,
                                             data=data2, tau.dist="HalfNormal",
                                             tau.prior=2.0, beta.prior=2,
                                             warmup=100, iter=200, chains=1, thin=1) )
              expect_true(nrow(fitted(map2)) == 4)
              data3 <- data.frame(n=c(n,n,n,n),r=c(5,5,5,5), study=c(1,1,2,2))
              suppressWarnings( map3 <- gMAP(cbind(r, n-r) ~ 1 | study, family=binomial,
                                             data=data3, tau.dist="HalfNormal",
                                             tau.prior=2.0, beta.prior=2,
                                             warmup=100, iter=200, chains=1, thin=1) )
              expect_true(nrow(fitted(map3)) == 4)
          })

test_that("gMAP handles fixed tau case", {
    suppressWarnings( map1 <- gMAP(cbind(r, n-r) ~ 1 | study, family=binomial,
                                   data=AS, tau.dist="Fixed",
                                   tau.prior=0.5, beta.prior=2,
                                   warmup=100, iter=200, chains=1, thin=1) )
    expect_true(map1$Rhat.max >= 1)
          })

test_that("gMAP labels data rows correctly when using covariates", {

    data_covs <- data.frame(n=10, r=3, study=c(1,2,2), stratum=factor(c("A", "A", "B")) ) %>%
        mutate(group=paste(study,stratum,sep="/"), id=as.integer(factor(group)))

    suppressMessages(suppressWarnings(map_covs <- gMAP(cbind(r, n-r) ~ 1 + stratum | study, family=binomial,
                                                       data=data_covs, tau.dist="Fixed",
                                                       tau.prior=0.25, beta.prior=2,
                                                       warmup=100, iter=200, chains=1, thin=1)))

    expect_true(all(rownames(fitted(map_covs)) == paste(data_covs$study, data_covs$stratum, sep="/")))

    suppressMessages(suppressWarnings(map_tau_strata <- gMAP(cbind(r, n-r) ~ 1 | id, family=binomial,
                                                             tau.strata=stratum,
                                                             data=data_covs, tau.dist="Fixed",
                                                             tau.prior=c(0.25, 0.5), beta.prior=2,
                                                             warmup=100, iter=200, chains=1, thin=1)))
    expect_true(all(rownames(fitted(map_tau_strata)) == as.character(data_covs$id)))

    })
