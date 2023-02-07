context("EM: Expectation-Maximization")

## test that the EM algorithms recover reliably test distributions;
## test criterium is a "sufficiently" small KL divergence

## the test considers the uni-variate normal, beta & gamma case in
## three flavours each:

## - single-component
## - two component mixture with heavy tails
## - three component mixture with bi-modal density and heavy tails


## number of samples drawn from test distributions
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    ## through testing if not on CRAN
    Nsim <- 1e4
    verbose <- FALSE
    KLthresh <- 1e-2
} else {
    ## on CRAN we shortcut
    Nsim <- 1e3
    verbose <- FALSE
    KLthresh <- 1e-1
}


## setup test cases

ref <- list()

ref$norm_single <- mixnorm(c(1, 1, 5),
                           param="mn", sigma=1)

ref$norm_heavy <- mixnorm(c(0.5, 0, 0.25),
                          c(0.5, 1, 5),
                          param="mn", sigma=1)

ref$norm_bi <- mixnorm(c(0.5, 0, 0.5),
                       c(0.25, 1, 5),
                       c(0.25, -1, 2),
                       param="mn", sigma=1)

ref$beta_single <- mixbeta(c(1, 0.3, 10),
                           param="mn")

## density which is challenging for the constrained version of the
## beta EM (and leads to a large KLdiv)
ref$beta_single_alt  <- mixbeta(c(1, 0.2, 3))

ref$beta_heavy <- mixbeta(c(0.8, 0.3, 10),
                            c(0.2, 0.5, 2.5),
                            param="mn")

ref$beta_bi <- mixbeta(c(0.3, 0.3, 20),
                         c(0.2, 0.5, 2),
                         c(0.5, 0.7, 10),
                         param="mn")

ref$gamma_single <- mixgamma(c(1, 7.5, 5),
                             param="mn",
                             likelihood="poisson")

ref$gamma_heavy <- mixgamma(c(0.5, 7.5, 0.5),
                            c(0.5, 5, 10),
                            param="mn",
                            likelihood="poisson")

ref$gamma_bi <- mixgamma(c(0.5, 7.5, 1),
                         c(0.25, 15, 15),
                         c(0.25, 5, 10),
                         param="mn",
                         likelihood="poisson")

EM_test <- function(mixTest, seed, Nsim=1e4, verbose=FALSE, ...) {
    set.seed(seed)
    samp <- rmix(mixTest, Nsim)
    EMmix <- mixfit(samp,
                    type=switch(class(mixTest)[1], gammaMix="gamma", normMix="norm", betaMix="beta"),
                    thin=1,
                    eps=2,
                    Nc=ncol(mixTest),
                    verbose=verbose, ...)
    kl <- abs(KLdivmix(mixTest, EMmix))
    expect_true(kl < KLthresh)
}

test_that("Normal EM fits single component",     EM_test(ref$norm_single, 3453563, Nsim, verbose))
test_that("Normal EM fits heavy-tailed mixture", EM_test(ref$norm_heavy,  9275624, Nsim, verbose))
test_that("Normal EM fits bi-modal mixture",     EM_test(ref$norm_bi,     9345726, Nsim, verbose))

test_that("Gamma EM fits single component",      EM_test(ref$gamma_single, 9345835, Nsim, verbose))
test_that("Gamma EM fits heavy-tailed mixture",  EM_test(ref$gamma_heavy,  5629389, Nsim, verbose))
test_that("Gamma EM fits bi-modal mixture",      EM_test(ref$gamma_bi,     9373515, Nsim, verbose))

test_that("Beta EM fits single component",       EM_test(ref$beta_single, 7265355, Nsim, verbose))
test_that("Beta EM fits single component with mass at boundary", EM_test(ref$beta_single_alt, 7265355, Nsim, verbose, constrain_gt1=FALSE))
test_that("Beta EM fits heavy-tailed mixture",   EM_test(ref$beta_heavy,  2946562, Nsim, verbose))
test_that("Beta EM fits bi-modal mixture",       EM_test(ref$beta_bi,     9460370, Nsim, verbose))

test_that("Constrained Beta EM respects a>1 & b>1", {
    unconstrained  <- mixbeta(c(0.6, 2.8, 64), c(0.25, 0.5, 0.92), c(0.15, 3, 15))
    set.seed(45747)
    samp <- rmix(unconstrained, Nsim)
    constrained  <- mixfit(samp, type="beta", Nc=3, constrain_gt1=TRUE)
    expect_numeric(constrained[2,], lower=1, any.missing=FALSE, len=3)
    expect_numeric(constrained[3,], lower=1, any.missing=FALSE, len=3)
}
)
