## One half of the RBesT.MC.s2z on/off pair.
##
## The switch changes the sampling geometry, not the model, so the two
## posteriors must agree up to Monte Carlo error. Both halves use the same
## seed, data, priors and sampler settings and differ only in RBesT.MC.s2z.
fixture <- withr::with_options(
  list(
    RBesT.MC.s2z = TRUE,
    RBesT.MC.save_warmup = FALSE
  ),
  {
    set.seed(5711)
    suppressMessages(suppressWarnings(
      gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = AS,
        family = binomial,
        tau.dist = "HalfNormal",
        tau.prior = 0.5,
        beta.prior = 2,
        warmup = 1000,
        iter = 2000,
        chains = 4,
        thin = 1
      )
    ))
  }
)
