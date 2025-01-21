## run the example from predict.gMAP
source_example("predict_gMAP.R")

## check that we got for each input data item a prediction
test_that("correct # of predictions are generated", {
  expect_equal(nrow(map$data), ncol(samp))
})

## check that the predictive distribution has a variance which is
## larger in accordance to the betwee-trial heterogeniety (needs to be
## done on the link scale)

test_that("variances have correct ordering", {
  pred_cov_link <- predict(map, type = "link")
  within_var <- (summary(pred_cov_link)[, "sd"])^2

  pred_cov_link_pred <- predict(map, trans_cov, type = "link")
  pred_var_pred <- summary(pred_cov_link_pred)[, "sd"]
  tau_est <- summary(map)$tau[, "mean"]

  ## the predictive must include between and within; as such it is
  ## larger than within
  expect_true(all(pred_var_pred > tau_est))

  ## ensure that predictive has larger variance than the model estimate
  expect_true(all(summary(pred_cov_link_pred)[, "sd"] > summary(pred_cov_link)[, "sd"]))
})


## new prediction was done for a single data item
test_that("correct # of new predictions are generated", {
  expect_equal(ncol(pred_new), 1)
})

## must have larger sd than between-trial alone (on link scale)
test_that("predictive variances have correct ordering", {
  pred_new_link <- predict(map, data.frame(country = "CH", study = 11), type = "link")
  tau_est <- summary(map)$tau[, "mean"]
  expect_true(summary(pred_new_link)[, "sd"] > tau_est)
})

## whenever the same study/covariate combination is requested, then
## the MAP must be numerically exactly the same. This ensures that per
## study the random effect is sampled just once in each iteration.
test_that("predictive distributions for the same study & covariate must match exactly", {
  trans_cov_new <- data.frame(study = "new", n = 50, r = 0, country = levels(trans_cov$country)[c(1, 1)])
  post_trans <- as.matrix(predict(map, newdata = trans_cov_new))
  expect_equal(post_trans[, 1], post_trans[, 2])
})

test_that("automixfit attempts K=4 different models and returns best fitting", {
  auto_map <- automixfit(map, Nc = 1:4, k = 6)
  models <- attr(auto_map, "models")
  expect_equal(length(models), 4)
  perf <- sapply(models, AIC, k = 6)
  ## ensure that performance is decreasing
  expect_true(all(diff(perf) > 0))
  expect_true("betaMix" %in% class(auto_map))
})


test_that("mixfit for prediction handles response and link scale", {
  pred_map <- mixfit(pred_new, Nc = 2)

  expect_true(is.list(pred_map))
  expect_true("betaMix" %in% class(pred_map[[1]]))
  expect_equal(ncol(pred_map[[1]]), 2)

  pred_new_link <- predict(map, data.frame(country = "CH", study = 11), type = "link")
  pred_map_link <- mixfit(pred_new_link, Nc = 2)

  expect_true(is.list(pred_map_link))
  expect_true("normMix" %in% class(pred_map_link[[1]]))
  expect_equal(ncol(pred_map_link[[1]]), 2)
})


source_example("mixcombine.R")

test_that("combination of mixtures", {
  m1 <- mixcombine(bm, unif, weight = c(9, 1))
  m2 <- mixcombine(bm, unif, unif, weight = c(8, 1, 1))
  expect_equal(m1[1, ], c(bm[1, ] - 0.1 / 2, 0.1), ignore_attr = TRUE)
  expect_equal(m1[2:3, 1:2], bm[2:3, 1:2], ignore_attr = TRUE)
  expect_equal(m2[2:3, 1:2], bm[2:3, 1:2], ignore_attr = TRUE)
})

test_that("throws an error if more weights than mixtures given", {
  ## giving 3 weights but only 2 mixtures must not work
  expect_error(mixcombine(bm, unif, weight = c(8, 1, 1)), "length(weight) not equal to length(comp)", fixed = TRUE)
})

test_that("combination of normal mixtures without default sigma works", {
  norm_ui <- mixnorm(c(1, 0, 2))
  norm_ui_mix <- mixcombine(norm_ui, norm_ui, weight = c(0.5, 0.5))
  expect_true(ncol(norm_ui_mix) == 2)
})

source_example("robustify.R")

test_that("beta mixture is robustified with Beta(1,1)", {
  expect_equal(ncol(bmix) + 1, ncol(rbmix))
  expect_equal(rbmix[, ncol(rbmix)], c(0.1, 1, 1), ignore_attr = TRUE)
})

test_that("beta mixture is robustified with Beta(0.5,0.5)", {
  rbmix2 <- robustify(bmix, w = 0.1, n = 0, mean = 0.5)
  expect_equal(ncol(bmix) + 1, ncol(rbmix2))
  expect_equal(rbmix2[, ncol(rbmix2)], c(0.1, 0.5, 0.5), ignore_attr = TRUE)
})

test_that("gamma mixture is robustified with n=1 equivalent prior", {
  m <- summary(gmnMix)["mean"]
  nr <- ncol(rgmnMix)
  expect_equal(rgmnMix[[nr, rescale = TRUE]], mixgamma(c(1, m, 1), param = "mn"), ignore_attr = TRUE)
  expect_equal(rgmnMix[1, nr], 0.1)
})

test_that("gamma mixture is robustified with n=5 equivalent prior", {
  m <- summary(gmnMix)["mean"]
  rgmnMix2 <- robustify(gmnMix, w = 0.1, n = 5, mean = 2)
  nr <- ncol(rgmnMix2)
  expect_equal(rgmnMix2[[nr, rescale = TRUE]], mixgamma(c(1, m, 5), param = "mn"), ignore_attr = TRUE)
  expect_equal(rgmnMix2[1, nr], 0.1)
})

test_that("normal mixture is robustified with n=1 equivalent prior", {
  nr <- ncol(rnMix)
  expect_equal(rnMix[[nr, rescale = TRUE]], mixnorm(c(1, 0, 1), param = "mn", sigma = sigma(nm)), ignore_attr = TRUE)
  expect_equal(rnMix[1, nr], 0.1)
})

test_that("normal mixture is robustified with n=5 equivalent prior", {
  rnMix2 <- robustify(nm, w = 0.1, mean = 0, n = 5, sigma = sigma(nm))
  nr <- ncol(rnMix2)
  expect_equal(rnMix2[[nr, rescale = TRUE]], mixnorm(c(1, 0, 5), param = "mn", sigma = sigma(nm)), ignore_attr = TRUE)
  expect_equal(rnMix2[1, nr], 0.1)
})

test_that("plotting of normal mixtures without default sigma works", {
  norm_ui <- mixnorm(c(1, 0, 2))
  norm_mix_ui <- mixcombine(norm_ui, norm_ui, weight = c(0.5, 0.5))
  pl <- plot(norm_mix_ui)
  expect_true(inherits(pl, "ggplot"))
})

source_example("ess.R")

test_that("conjugate beta case matches canonical formula", {
  expect_equal(a + b, ess(prior, "moment"))
  expect_equal(a + b, round(ess(prior, "morita")))
  expect_equal(a + b, ess(prior, "elir"))
})

test_that("ess elir for beta mixtures gives a warning for a<1 & b<1 densities", {
  unconstrain1 <- mixbeta(c(0.95, 10, 5), c(0.05, 0.9, 2))
  unconstrain2 <- mixbeta(c(0.95, 10, 5), c(0.05, 2, 0.9))

  expect_error(ess(unconstrain1, "elir"), "At least one parameter of the beta mixtures is less than 1")
  expect_error(ess(unconstrain2, "elir"), "At least one parameter of the beta mixtures is less than 1")

  ## this one can trigger errors if the integration is not setup properly
  constrained <- mixbeta(c(0.48, 1, 11), c(0.34, 6.9, 173), c(0.18, 1.0, 1.13))
  expect_numeric(ess(constrained, "elir"), lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
})

test_that("conjugate normal case matches canonical formula", {
  s <- 2
  sigma_data <- 4
  nprior <- mixnorm(c(1, -1, s), sigma = sigma_data)
  nprior_ess <- sigma_data^2/s^2
  expect_equal(ess(nprior, "moment", sigma = sigma_data), nprior_ess)
  expect_equal(ess(nprior, "morita", sigma = sigma_data, s=Inf), nprior_ess)
  expect_equal(ess(nprior, "elir", sigma = sigma_data), nprior_ess)
})

test_that("ess elir for normal mixtures returns correct values", {
  mix <- mixnorm(inf1 = c(0.5026, -191.1869, 127.4207), inf2 = c(0.2647, -187.5895, 31.6130), inf3 = c(0.2326, -184.7445, 345.3849), sigma = 270.4877)
  expect_gt(ess(mix, sigma = 270.4877), 0)
})

test_that("moment matching for beta mixtures is correct", {
  expect_equal(ess(bmix, method = "moment"), sum(ab_matched))
})

test_that("beta mix ess works when run through sapply", {
  expect_numeric(sapply(X = list(bmix), FUN=ess))
})

test_that("normal mixtures have reference scale used correctly", {
  nmix_sigma_small <- nmix
  nmix_sigma_large <- nmix
  sigma_large <- 2 * summary(nmix_sigma_large)["sd"]
  sigma(nmix_sigma_large) <- sigma_large
  sigma(nmix_sigma_small) <- sigma_large / sqrt(2)
  suppressMessages(e1m <- ess(nmix_sigma_large, "moment"))
  suppressMessages(e2m <- ess(nmix_sigma_small, "moment"))
  expect_gt(e1m, e2m)
  expect_equal(floor(abs(e2m - e1m / 2)), 0)

  suppressMessages(e1b <- ess(nmix_sigma_large, "morita"))
  suppressMessages(e2b <- ess(nmix_sigma_small, "morita"))
  expect_gt(e1b, e2b)
  expect_equal(floor(abs(e2b - e1b / 2)), 0)

  suppressMessages(e1r <- ess(nmix_sigma_large, "elir"))
  suppressMessages(e2r <- ess(nmix_sigma_small, "elir"))
  expect_gt(e1r, e2r)
  expect_equal(floor(abs(e2r - e1r / 2)), 0)
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

test_that("gamma mix ess works when run through sapply", {
  expect_numeric(sapply(X = list(gmix), FUN=ess))
})

test_that("gamma 1-component density gives canonical results", {
  guni1 <- gmix[[1, rescale = TRUE]]
  likelihood(guni1) <- "poisson"
  guni2 <- gmix[[1, rescale = TRUE]]
  likelihood(guni2) <- "exp"

  e1m <- ess(guni1, "moment")
  e2m <- ess(guni2, "moment")
  expect_true(e1m != e2m)
  expect_equal(guni1[3, 1], e1m)
  expect_equal(guni2[2, 1], e2m)

  e1b <- round(ess(guni1, "morita"))
  e2b <- round(ess(guni2, "morita"))
  expect_true(e1b != e2b)
  expect_equal(guni1[3, 1], e1b)
  expect_equal(guni2[2, 1], e2b)

  e1r <- ess(guni1, "elir")
  e2r <- ess(guni2, "elir")
  expect_true(e1r != e2r)
  expect_true(abs(guni1[3, 1] - e1r) < 1E-4)
  ## ELIR gives a-1 as ESS
  expect_true(abs(guni2[2, 1] - (e2r + 1)) < 1E-4)
})

## check predictive consistency of ELIR
elir_predictive_consistent <- function(dens, m, Nsim, seed, stat, ...) {
  ## simulated from predictve which is m events equivalent to
  suppressMessages(pdens <- preddist(dens, n = m))
  set.seed(seed)
  psamp <- rmix(pdens, Nsim)

  if (inherits(dens, "gammaMix")) {
    psamp <- psamp / m
  }

  posterior_ess <- function(mix, method, stat, ...) {
    args <- c(list(priormix = mix, stat = 0), list(...))
    names(args)[2] <- stat
    fn <- function(x) {
      args[[stat]] <- x
      suppressMessages(res <- ess(do.call(postmix, args), method = method))
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
  nmix <- mixnorm(rob = c(0.5, 0, 2), inf = c(0.5, 3, 4), sigma = 10)
  elir_predictive_consistent(nmix, m = 3E2, Nsim = 1E3, seed = 3435, stat = "m", se = 10 / sqrt(3E2))
})

test_that("ESS elir is predictively consistent for beta mixtures", {
  skip_on_cran()
  bmix <- mixbeta(rob = c(0.2, 1, 1), inf = c(0.8, 10, 2))
  elir_predictive_consistent(bmix, m = 1E2, Nsim = 1E3, seed = 355435, stat = "r", n = 1E2)
})

test_that("ESS elir is predictively consistent for gamma mixtures (Poisson likelihood)", {
  skip_on_cran()
  gmixP <- mixgamma(rob = c(0.3, 20, 4), inf = c(0.7, 50, 10), likelihood = "poisson")
  elir_predictive_consistent(gmixP, m = 1E2, Nsim = 1E3, seed = 355435, stat = "m", n = 1E2)
})

test_that("ess elir for problematic beta mixtures gives correct result 1", {
  ## by user reported beta mixture density which triggers this erros
  ## with RBesT 1.7.2 & 1.7.3 (others not tested):
  ## Error in if (all(dgl < 0) || all(dgl > 0)) { :
  ##   missing value where TRUE/FALSE needed

  mixmat <- matrix(c(
    0.06429517, 0.03301215, 0.00269268, 0.90000000,
    437.32302999, 64.04211307, 5.92543558, 1.00000000,
    10.71709277, 2.14157953, 1.00000001, 1.00000000
  ), byrow = TRUE, ncol = 4)

  mixb <- do.call(mixbeta, apply(mixmat, 2, c, simplify = FALSE))

  expect_double(ess(mixb), lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
})

test_that("ess elir for problematic beta mixtures gives correct result 2", {
  mixmat <- matrix(c(
    0.7237396, 0.1665037, 0.1097567,
    53.3721902, 44.3894573, 9.8097062,
    1.4301638, 4.3842200, 1.8492197
  ), byrow = TRUE, ncol = 3)

  mixb <- do.call(mixbeta, apply(mixmat, 2, c, simplify = FALSE))

  expect_double(ess(robustify(mixb, 0.05, 0.5)), lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  expect_double(ess(robustify(mixb, 0.95, 0.5)), lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
})


test_that("ess elir for problematic beta mixtures gives warning", {
  mixmat1 <- matrix(c(
    0.6092774, 0.2337629, 0.1569597,
    1.0000000, 1.2672179, 3.3856153,
    11.8465288, 1.2389927, 7.0191159
  ), byrow = TRUE, ncol = 3)


  mixb1 <- do.call(mixbeta, apply(mixmat1, 2, c, simplify = FALSE))

  ## in case one of the coefficients of a and b is 1, then we can
  ## get negative results... which are unreliable to the user hopefully
  expect_warning(ess(mixb1))
  expect_double(suppressWarnings(ess(mixb1)), finite = TRUE,
                any.missing = FALSE, len = 1)

  mixmat2 <- matrix(c(
    0.6051804, 0.2324492, 0.1623704,
    1.0210697, 1.1955047, 3.1342298,
    11.5485831, 1.0831573, 6.7636286
  ), byrow = TRUE, ncol = 3)

  mixb2 <- do.call(mixbeta, apply(mixmat2, 2, c, simplify = FALSE))

  expect_double(ess(mixb2), lower = 0, finite = TRUE, any.missing = FALSE,
                len = 1)
})

test_that("BinaryExactCI has correct boundary behavior", {
  expect_equal(unname(BinaryExactCI(0, 10, 0.05)[1]), 0)
  expect_equal(unname(BinaryExactCI(10, 10, 0.05)[2]), 1)
})


test_that("ess for a normal density with binomial family under a logit link gives correct results", {
  ## the ess elir for a normal density prior with mean m and
  ## standard deviation s given to a logit transformed response rate
  ## is: i(p(eta)) = 1/s^2 and i_F(eta) = exp(eta) / (1 +
  ## exp(eta))^2 => r(eta) = i(p(eta)) / i_F(eta) the ess ELIR
  ## integral then involves terms as integral exp(eta) p(eta|m,s)
  ## d(eta) = exp(m + s^2/2) and integral exp(-eta) p(eta|m,s) d(eta)
  ## = exp(-m + s^2/2). The analyical result is then
  ## ess_elir = 1/s^2 * [ 2 + exp(-m + s^/2) + exp(m + s^/2) ]

  ## since the information of a normal is just 1/s^2 and thus a
  ## constant, the moment based approach gives the same result as
  ## the elir method. The morita method differs though as it
  ## evaluates at the mode of the prior.

  ess_elir_binomial_logit <- function(m, s) {
    s2 <- s * s
    (2 + exp(-m + s2 / 2) + exp(m + s2 / 2)) / s2
  }

  fisher_binomial_logit <- function(l) {
    exp(l) / (1 + exp(l))^2
  }

  ## Pennello and Thomson ESS for a binomial logit case corresponds
  ## to the Morita ESS whenever the scale s for the flattened prior
  ## is Infinity
  pe_ess_binomial_logit <- function(m, s) {
    1 / s^2 / fisher_binomial_logit(m)
  }

  m1 <- 0
  s1 <- 2 / sqrt(10)
  prior_norm1 <- mixnorm(c(1, m1, s1), sigma = 2)
  expect_equal(ess(prior_norm1, "elir", family = binomial, sigma = 2), ess_elir_binomial_logit(m1, s1), tolerance = 1E-4)
  expect_equal(ess(prior_norm1, "moment", family = binomial, sigma = 2), ess_elir_binomial_logit(m1, s1), tolerance = 1E-4)
  expect_equal(ess(prior_norm1, "morita", family = binomial, sigma = 2, s = Inf), pe_ess_binomial_logit(m1, s1), tolerance = 1E-4)

  m2 <- 2
  s2 <- 2 / sqrt(100)
  prior_norm2a <- mixnorm(c(1, m2, s2), sigma = 2)
  expect_equal(ess(prior_norm2a, "elir", family = binomial, sigma = 2), ess_elir_binomial_logit(m2, s2), tolerance = 1E-4)
  expect_equal(ess(prior_norm2a, "moment", family = binomial, sigma = 2), ess_elir_binomial_logit(m2, s2), tolerance = 1E-4)
  expect_equal(ess(prior_norm2a, "morita", family = binomial, sigma = 2, s = Inf), pe_ess_binomial_logit(m2, s2), tolerance = 1E-4)

  ## sigma does not play a role here
  prior_norm2b <- mixnorm(c(1, m2, s2), sigma = 4)
  expect_equal(ess(prior_norm2b, "elir", family = binomial, sigma = 4), ess_elir_binomial_logit(m2, s2), tolerance = 1E-4)
  expect_equal(ess(prior_norm2b, "moment", family = binomial, sigma = 4), ess_elir_binomial_logit(m2, s2), tolerance = 1E-4)
  expect_equal(ess(prior_norm2b, "morita", family = binomial, sigma = 4, s = Inf), pe_ess_binomial_logit(m2, s2), tolerance = 1E-4)
})


test_that("ess for a normal density with poisson family under a log link gives correct results", {
  ## the ess elir for a normal density prior with mean m and
  ## standard deviation s given to a log transformed count rate
  ## is: i(p(eta)) = 1/s^2 and i_F(eta) = exp(eta)  => r(eta) = i(p(eta)) / i_F(eta) the ess ELIR
  ## integral then involves terms as integral exp(-eta) p(eta|m,s)
  ## d(eta) = exp(-m + s^2/2). The analyical result is then
  ## ess_elir = 1/s^2 * exp(-m + s^/2)

  ## since the information of a normal is just 1/s^2 and thus a
  ## constant, the moment based approach gives the same result as
  ## the elir method. The morita method differs though as it
  ## evaluates at the mode of the prior.

  ess_elir_poisson_log <- function(m, s) {
    s2 <- s * s
    exp(-m + s2 / 2) / s2
  }

  fisher_poisson_log <- function(l) {
    exp(l)
  }

  ## Pennello and Thomson ESS for a binomial logit case corresponds
  ## to the Morita ESS whenever the scale s for the flattened prior
  ## is Infinity
  pe_ess_poisson_log <- function(m, s) {
    1 / s^2 / fisher_poisson_log(m)
  }

  m1 <- 0
  s1 <- 2 / sqrt(10)
  prior_norm1 <- mixnorm(c(1, m1, s1), sigma = 2)
  expect_equal(ess(prior_norm1, "elir", family = poisson, sigma = 2), ess_elir_poisson_log(m1, s1), tolerance = 1E-4)
  expect_equal(ess(prior_norm1, "moment", family = poisson, sigma = 2), ess_elir_poisson_log(m1, s1), tolerance = 1E-4)
  expect_equal(ess(prior_norm1, "morita", family = poisson, sigma = 2, s = Inf), pe_ess_poisson_log(m1, s1), tolerance = 1E-4)

  m2 <- 2
  s2 <- 2 / sqrt(100)
  prior_norm2a <- mixnorm(c(1, m2, s2), sigma = 2)
  expect_equal(ess(prior_norm2a, "elir", family = poisson, sigma = 2), ess_elir_poisson_log(m2, s2), tolerance = 1E-4)
  expect_equal(ess(prior_norm2a, "moment", family = poisson, sigma = 2), ess_elir_poisson_log(m2, s2), tolerance = 1E-4)
  expect_equal(ess(prior_norm2a, "morita", family = poisson, sigma = 2, s = Inf), pe_ess_poisson_log(m2, s2), tolerance = 1E-4)

  ## sigma does not play a role here
  prior_norm2b <- mixnorm(c(1, m2, s2), sigma = 4)
  expect_equal(ess(prior_norm2b, "elir", family = poisson, sigma = 4), ess_elir_poisson_log(m2, s2), tolerance = 1E-4)
  expect_equal(ess(prior_norm2b, "moment", family = poisson, sigma = 4), ess_elir_poisson_log(m2, s2), tolerance = 1E-4)
  expect_equal(ess(prior_norm2b, "morita", family = poisson, sigma = 4, s = Inf), pe_ess_poisson_log(m2, s2), tolerance = 1E-4)
})

test_that("ess for a beta mixture errors if a family is specified", {
  expect_error(ess(mixbeta(c(1, 5, 15)), family = binomial))
})

test_that("ess for a gamma mixture errors if a family is specified", {
  expect_error(ess(mixgamma(rob = c(0.3, 20, 4), inf = c(0.7, 50, 10)), family = poisson))
})
