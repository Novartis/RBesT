context("postmix: Posterior Mixture Distribution")

norm <- mixnorm(c(1, 0, 0.5), sigma=1)

test_that("Normal mixture reference scale is updated", {
              suppressMessages(post_norm <- postmix(norm, m=0, n=100, se=0.1))
              expect_equal(1, RBesT::sigma(post_norm))
          })

test_that("Normal mixture default reference scale is used", {
              suppressMessages(post_norm <- postmix(norm, m=0, n=100))
              psd <- sqrt( 1/( (1/sqrt(100))^-2 + (1/2)^-2) )
              expect_lt(abs(summary(post_norm)["sd"] - psd), 1E-7)
          })

test_that("Normal mixture default reference scale is updated", {
              suppressMessages(post_norm <- postmix(norm, m=0, se=1, n=100))
              expect_equal(sigma(post_norm), 10)
          })


test_that("Gamma mixture is updated for Poisson likelihood", {
              gamma_prior <- mixgamma(c(1,10,1), param="mn")
              gamma_post  <- postmix(gamma_prior, n=20, m=2)
              expect_equal(ess(gamma_post), 21)
          })
