context("Utilities: predict")

## run the example from predict.gMAP
suppressMessages(example("predict.gMAP", package="RBesT", echo=FALSE, ask=FALSE, verbose=FALSE))

## check that we got for each input data item a prediction
test_that("correct # of predictions are generated", expect_equal(nrow(map$data), ncol(samp)))

## check that the predictive distribution has a variance which is
## larger in accordance to the betwee-trial heterogeniety (needs to be
## done on the link scale)

test_that("variances have correct ordering", {
              pred_cov_link <- predict(map, type="link")
              within_var <- (summary(pred_cov_link)[,"sd"])^2

              pred_cov_link_pred <- predict(map, trans_cov, type="link")
              pred_var_pred <- summary(pred_cov_link_pred)[,"sd"]
              tau_est <- summary(map)$tau[,"mean"]

              ## the predictive must include between and within; as such it is
              ## larger than within
              expect_true(all(pred_var_pred > tau_est))

              ## ensure that predictive has larger variance than the model estimate
              expect_true(all(summary(pred_cov_link_pred)[,"sd"] > summary(pred_cov_link)[,"sd"]))
          })


## new prediction was done for a single data item
test_that("correct # of new predictions are generated", expect_equal(ncol(pred_new), 1))

## must have larger sd than between-trial alone (on link scale)
test_that("predictive variances have correct ordering",{
              pred_new_link <- predict(map, data.frame(country="CH", study=11), type="link")
              tau_est <- summary(map)$tau[,"mean"]
              expect_true(summary(pred_new_link)[,"sd"] > tau_est)
          })

## whenever the same study/covariate combination is requested, then
## the MAP must be numerically exactly the same. This ensures that per
## study the random effect is sampled just once in each iteration.
test_that("predictive distributions for the same study & covariate must match exactly", {
    trans_cov_new <- data.frame(study="new", n=50, r=0, country=levels(trans_cov$country)[c(1,1)])
    post_trans <- as.matrix(predict(map, newdata=trans_cov_new))
    expect_equal(post_trans[,1], post_trans[,2])
})

context("Utilities: (auto)mixfit")

test_that("automixfit attempts K=4 different models and returns best fitting", {
              auto_map <- automixfit(map, Nc=1:4, k=6)
              models <- attr(auto_map, "models")
              expect_equal(length(models), 4)
              perf <- sapply(models, AIC, k=6)
              ## ensure that performance is decreasing
              expect_true(all(diff(perf) > 0))
              expect_true("betaMix" %in% class(auto_map))
          })


test_that("mixfit for prediction handles response and link scale", {
              pred_map <- mixfit(pred_new, Nc=2)

              expect_true(is.list(pred_map))
              expect_true("betaMix" %in% class(pred_map[[1]]))
              expect_equal(ncol(pred_map[[1]]), 2)

              pred_new_link <- predict(map, data.frame(country="CH", study=11), type="link")
              pred_map_link <- mixfit(pred_new_link, Nc=2)

              expect_true(is.list(pred_map_link))
              expect_true("normMix" %in% class(pred_map_link[[1]]))
              expect_equal(ncol(pred_map_link[[1]]), 2)
          })


context("Utilities: mixcombine")

example("mixcombine", package="RBesT", echo=FALSE, ask=FALSE, verbose=FALSE)

test_that("combination of mixtures", {
              m1 <- mixcombine(bm, unif, weight=c(9, 1))
              m2 <- mixcombine(bm, unif, unif, weight=c(8, 1, 1))
              expect_equivalent(m1[1,], c(bm[1,] - 0.1/2, 0.1))
              expect_equivalent(m1[2:3,1:2], bm[2:3,1:2])
              expect_equivalent(m2[2:3,1:2], bm[2:3,1:2])
          })

test_that("throws an error if more weights than mixtures given", {
              ## giving 3 weights but only 2 mixtures must not work
              expect_error(mixcombine(bm, unif, weight=c(8, 1, 1)), "length(weight) not equal to length(comp)", fixed=TRUE)
          })

test_that("combination of normal mixtures without default sigma works", {
              norm_ui <- mixnorm(c(1, 0, 2))
              norm_ui_mix <- mixcombine(norm_ui, norm_ui, weight=c(0.5,0.5))
              expect_true(ncol(norm_ui_mix) == 2)
          })

context("Utilities: robustify")

example("robustify", package="RBesT", echo=FALSE, ask=FALSE, verbose=FALSE)

test_that("beta mixture is robustified with Beta(1,1)", {
              expect_equal(ncol(bmix)+1, ncol(rbmix))
              expect_equivalent(rbmix[,ncol(rbmix)], c(0.1, 1, 1))
          })

test_that("beta mixture is robustified with Beta(0.5,0.5)", {
              rbmix2 <- robustify(bmix, w=0.1, n=0)
              expect_equal(ncol(bmix)+1, ncol(rbmix2))
              expect_equivalent(rbmix2[,ncol(rbmix2)], c(0.1, 0.5, 0.5))
          })

test_that("gamma mixture is robustified with n=1 equivalent prior", {
              m <- summary(gmnMix)["mean"]
              nr <- ncol(rgmnMix)
              expect_equivalent(rgmnMix[[nr, rescale=TRUE]], mixgamma(c(1, m, 1), param="mn"))
              expect_equal(rgmnMix[1,nr], 0.1)
          })

test_that("gamma mixture is robustified with n=5 equivalent prior", {
              m <- summary(gmnMix)["mean"]
              rgmnMix2 <- robustify(gmnMix, w=0.1, n=5)
              nr <- ncol(rgmnMix2)
              expect_equivalent(rgmnMix2[[nr, rescale=TRUE]], mixgamma(c(1, m, 5), param="mn"))
              expect_equal(rgmnMix2[1,nr], 0.1)
          })

test_that("normal mixture is robustified with n=1 equivalent prior", {
              nr <- ncol(rnMix)
              expect_equivalent(rnMix[[nr, rescale=TRUE]], mixnorm(c(1, 0, 1), param="mn", sigma=sigma(nm)))
              expect_equal(rnMix[1,nr], 0.1)
          })

test_that("normal mixture is robustified with n=5 equivalent prior", {
              rnMix2 <- robustify(nm, w=0.1, mean=0, n=5, sigma=sigma(nm))
              nr <- ncol(rnMix2)
              expect_equivalent(rnMix2[[nr, rescale=TRUE]], mixnorm(c(1, 0, 5), param="mn", sigma=sigma(nm)))
              expect_equal(rnMix2[1,nr], 0.1)
          })

context("Utilities: Plotting of Mixtures")
test_that("plotting of normal mixtures without default sigma works", {
              norm_ui <- mixnorm(c(1, 0, 2))
              norm_mix_ui <- mixcombine(norm_ui, norm_ui, weight=c(0.5,0.5))
              pl <- plot(norm_mix_ui)
              expect_true(inherits(pl, "ggplot"))
          })

context("Utilities: Mixture Effective Sample Size")

example("ess", package="RBesT", echo=FALSE, ask=FALSE, verbose=FALSE)

test_that("conjugate beta case matches canonical formula", {
              expect_equal(a+b, ess(prior, "moment"))
              expect_equal(a+b, round(ess(prior, "morita")))
              expect_equal(a+b, ess(prior, "elir"))
          })

test_that("ess elir for beta mixtures gives a warning for a<1 & b<1 densities", {
    unconstrain1 <- mixbeta(c(0.95, 10, 5), c(0.05, 0.9, 2))
    unconstrain2 <- mixbeta(c(0.95, 10, 5), c(0.05, 2, 0.9))

    expect_error(ess(unconstrain1, "elir"), "At least one parameter of the beta mixtures is less than 1")
    expect_error(ess(unconstrain2, "elir"), "At least one parameter of the beta mixtures is less than 1")

    ## this one can trigger errors if the integration is not setup properly
    constrained  <- mixbeta(c(0.48, 1, 11), c(0.34, 6.9, 173), c(0.18, 1.0, 1.13))
    expect_numeric(ess(constrained, "elir"), lower=0, finite=TRUE, any.missing=FALSE, len=1)
})

test_that("ess elir for normal mixtures returns correct values", {
    mix <- mixnorm( inf1=c(0.5026,-191.1869,127.4207),inf2=c(0.2647,-187.5895,31.6130),inf3=c(0.2326,-184.7445,345.3849), sigma=270.4877)
    expect_gt(ess(mix), 0)
})

test_that("moment matching for beta mixtures is correct", {
              expect_equal(ess(bmix, method="moment"), sum(ab_matched))
          })

test_that("normal mixtures have reference scale used correctly", {
              nmix_sigma_small <- nmix
              sigma_large <- RBesT::sigma(nmix)
              sigma(nmix_sigma_small) <- sigma_large/sqrt(2)
              suppressMessages(e1m <- ess(nmix, "moment"))
              suppressMessages(e2m <- ess(nmix_sigma_small, "moment"))
              expect_gt(e1m, e2m)
              expect_equal(floor(abs(e2m - e1m/2)), 0)

              suppressMessages(e1b <- ess(nmix, "morita"))
              suppressMessages(e2b <- ess(nmix_sigma_small, "morita"))
              expect_gt(e1b, e2b)
              expect_equal(floor(abs(e2b - e1b/2)), 0)

              suppressMessages(e1r <- ess(nmix, "elir"))
              suppressMessages(e2r <- ess(nmix_sigma_small, "elir"))
              expect_gt(e1r, e2r)
              expect_equal(floor(abs(e2r - e1r/2)), 0)
          })

test_that("gamma mixtures have likelihood property respected", {
              gmix1 <- gmix
              likelihood(gmix1) <- "poisson"
              gmix2 <- gmix
              likelihood(gmix2) <- "exp"
              e1m <- ess(gmix1, "moment")
              e2m <- ess(gmix2, "moment")
              expect_true(e1m != e2m)

              e1b <- ess(gmix1, "morita")
              e2b <- ess(gmix2, "morita")
              expect_true(e1b != e2b)

              e1r <- ess(gmix1, "morita")
              e2r <- ess(gmix2, "morita")
              expect_true(e1r != e2r)
          })


test_that("gamma 1-component density gives canonical results", {
              guni1 <- gmix[[1, rescale=TRUE]]
              likelihood(guni1) <- "poisson"
              guni2 <- gmix[[1, rescale=TRUE]]
              likelihood(guni2) <- "exp"

              e1m <- ess(guni1, "moment")
              e2m <- ess(guni2, "moment")
              expect_true(e1m != e2m)
              expect_equal(guni1[3,1], e1m)
              expect_equal(guni2[2,1], e2m)

              e1b <- round(ess(guni1, "morita"))
              e2b <- round(ess(guni2, "morita"))
              expect_true(e1b != e2b)
              expect_equal(guni1[3,1], e1b)
              expect_equal(guni2[2,1], e2b)

              e1r <- ess(guni1, "elir")
              e2r <- ess(guni2, "elir")
              expect_true(e1r != e2r)
              expect_true(abs(guni1[3,1] - e1r) < 1E-4)
              ## ELIR gives a-1 as ESS
              expect_true(abs(guni2[2,1] - (e2r+1)) < 1E-4)
          })

## check predictive consistency of ELIR
elir_predictive_consistent  <- function(dens, m, Nsim, seed, stat, ...) {
    ## simulated from predictve which is m events equivalent to
    suppressMessages(pdens <- preddist(dens, n=m))
    set.seed(seed)
    psamp  <- rmix(pdens, Nsim)

    if(inherits(dens, "gammaMix"))
        psamp <- psamp / m

    posterior_ess  <- function(mix, method, stat, ...) {
        args <- c(list(priormix=mix, stat=0), list(...))
        names(args)[2] <- stat
        fn  <- function(x) {
            args[[stat]] <- x
            suppressMessages(res  <- ess(do.call(postmix, args), method=method))
            res
        }
        Vectorize(fn)
    }

    ## obtain ess of each posterior
    pred_ess <- posterior_ess(dens, "elir", stat, ...)
    ess_psamp <- pred_ess(psamp)

    suppressMessages(elir_prior <- ess(dens, "elir"))
    ## the average over the predicitve of the posterior ESS must match
    ## the the elir value taken directly (when m is subtracted, of
    ## course)
    elir_pred <- mean(ess_psamp) - m

    expect_true(abs(elir_prior - elir_pred) < 0.75)
}


test_that("ESS elir is predictively consistent for normal mixtures", {
    skip_on_cran()
    nmix <- mixnorm(rob=c(0.5, 0, 2), inf=c(0.5, 3, 4), sigma=10)
    elir_predictive_consistent(nmix, m=3E2, Nsim=1E3, seed=3435, stat="m", se=10/sqrt(3E2))
})

test_that("ESS elir is predictively consistent for beta mixtures", {
    skip_on_cran()
    bmix <- mixbeta(rob=c(0.2, 1, 1), inf=c(0.8, 10, 2))
    elir_predictive_consistent(bmix, m=1E2, Nsim=1E3, seed=355435, stat="r", n=1E2)
})

test_that("ESS elir is predictively consistent for gamma mixtures (Poisson likelihood)", {
    skip_on_cran()
    gmixP <- mixgamma(rob=c(0.3, 20, 4), inf=c(0.7, 50, 10), likelihood="poisson")
    elir_predictive_consistent(gmixP, m=1E2, Nsim=1E3, seed=355435, stat="m", n=1E2)
})
