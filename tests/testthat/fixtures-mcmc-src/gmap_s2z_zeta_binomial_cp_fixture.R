## Verbose S2Z fit used by the zeta invariant under the centered parametrization.
##
## `RBesT.verbose = TRUE` keeps the raw parameters (`beta_raw`, `tau_raw`,
## `xi_eta`, `xi_abar`) in the draws, which is what makes `xi_abar[1]`
## observable at all.
fixture <- withr::with_options(
  list(
    RBesT.verbose = TRUE,
    RBesT.MC.ncp = 0,
    RBesT.MC.save_warmup = FALSE
  ),
  {
    set.seed(46711)
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
