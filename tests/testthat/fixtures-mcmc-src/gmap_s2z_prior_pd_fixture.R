## Prior predictive S2Z fit: the sharpest assertion in the S2Z test suite.
##
## Under `prior_PD` the recovered intercept marginal is known in closed form
## whatever tau does, `beta[1] ~ Normal(m1, s1)`, because the widening of the
## alpha prior and the recovery of abar must compose back exactly. Splitting
## the recovery as r instead of r^2, or mismatching the widened variance,
## breaks the sd while leaving every structural invariant intact.
##
## `tau.prior = 2` is chosen so that `tau/sqrt(J)` is comparable to `s1`: the
## recovered share `r = sd_a/sd_alpha` is then well away from 0, which is what
## gives the test power. With a tight tau prior r is tiny, abar is nearly zero
## and almost any recovery formula would pass.
##
## The keys `m1 = 0.7` and `s1 = 0.5` are re-derived in the test from
## `fixture$fit.data$beta_prior`, so they cannot drift apart from this recipe.
fixture <- withr::with_options(
  list(RBesT.MC.save_warmup = FALSE),
  {
    set.seed(9134)
    suppressMessages(suppressWarnings(
      gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = AS,
        family = binomial,
        tau.dist = "HalfNormal",
        tau.prior = 2,
        beta.prior = cbind(0.7, 0.5),
        prior_PD = TRUE,
        warmup = 1000,
        iter = 5000,
        chains = 4,
        thin = 1
      )
    ))
  }
)
