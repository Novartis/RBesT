# test S3 methods in alphabetical order
test_that("as_draws and friends have resonable outputs", {
  skip_on_cran()

  n_iter <- 200
  n_warmup <- 100
  n_chains <- 2
  suppressMessages(suppressWarnings(
    {
      withr::with_options(list(RBesT.MC.save_warmup = FALSE), {
        set.seed(34563)
        map <- gMAP(
          cbind(r, n - r) ~ 1 | study,
          family = binomial,
          data = AS,
          tau.dist = "Fixed",
          tau.prior = 0.5,
          beta.prior = 2,
          warmup = n_warmup,
          iter = n_iter,
          chains = n_chains,
          thin = 1
        )
      })
      withr::with_options(list(RBesT.MC.save_warmup = TRUE), {
        set.seed(34563)
        map_full <- gMAP(
          cbind(r, n - r) ~ 1 | study,
          family = binomial,
          data = AS,
          tau.dist = "Fixed",
          tau.prior = 0.5,
          beta.prior = 2,
          warmup = n_warmup,
          iter = n_iter,
          chains = n_chains,
          thin = 1
        )
      })
    }
  ))

  draws <- as_draws(
    map,
    variable = "theta_resp_pred"
  )
  expect_s3_class(draws, "draws_list")
  expect_equal(
    posterior::variables(draws),
    "theta_resp_pred"
  )
  expect_equal(posterior::ndraws(draws), nsamples(map))

  draws <- suppressMessages(as_draws_matrix(
    map,
    variable = "theta_resp_pred"
  ))
  expect_s3_class(draws, "draws_matrix")
  expect_equal(
    posterior::variables(draws),
    "theta_resp_pred"
  )
  expect_equal(posterior::ndraws(draws), nsamples(map))

  draws <- as_draws_array(
    map,
    variable = "theta_resp_pred"
  )
  expect_s3_class(draws, "draws_array")
  expect_equal(
    posterior::variables(draws),
    "theta_resp_pred"
  )
  expect_equal(posterior::ndraws(draws), nsamples(map))

  draws <- as_draws_df(
    map,
    variable = "theta_resp_pred"
  )
  expect_s3_class(draws, "draws_df")
  expect_equal(
    posterior::variables(draws),
    "theta_resp_pred"
  )
  expect_equal(posterior::ndraws(draws), nsamples(map))

  draws <- as_draws_list(
    map,
    variable = "theta_resp_pred"
  )
  expect_s3_class(draws, "draws_list")
  expect_equal(
    posterior::variables(draws),
    "theta_resp_pred"
  )
  expect_equal(posterior::ndraws(draws), nsamples(map))

  draws <- as_draws_rvars(map)
  expect_s3_class(draws, "draws_rvars")
  expect_true(posterior::nvariables(draws) > 0)
  expect_equal(posterior::ndraws(draws), nsamples(map))
  n_saved_samples <- sum(map$fit@sim$n_save)

  expect_equal(posterior::ndraws(draws), n_saved_samples)
  expect_equal(posterior::ndraws(draws), (n_iter - n_warmup) * n_chains)
  expect_equal(nsamples(map), (n_iter - n_warmup) * n_chains)

  draws_full <- as_draws_rvars(map_full, inc_warmup = TRUE)
  expect_s3_class(draws_full, "draws_rvars")
  expect_true(posterior::nvariables(draws_full) > 0)
  expect_equal(
    posterior::ndraws(draws_full),
    nsamples(map_full) + n_warmup * n_chains
  )
  n_saved_samples_full <- sum(map_full$fit@sim$n_save)

  expect_equal(posterior::ndraws(draws_full), n_saved_samples_full)
  expect_equal(posterior::ndraws(draws_full), n_iter * n_chains)
  expect_equal(nsamples(map_full), (n_iter - n_warmup) * n_chains)
})
