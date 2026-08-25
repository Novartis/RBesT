## One half of the adapt_delta wiring pair.
##
## `.gmap_sampler_control()` is unit tested directly, but that only proves the
## list is built correctly. This pair proves the value actually reaches the
## sampler: a higher target acceptance rate must produce a smaller adapted step
## size. The two halves are identical apart from adapt_delta, and both hold the
## parametrization fixed, since step sizes are only comparable within one
## geometry.
fixture <- withr::with_options(
  list(
    RBesT.MC.ncp = 1,
    RBesT.MC.s2z = TRUE,
    RBesT.MC.control = list(adapt_delta = 0.99),
    RBesT.MC.save_warmup = FALSE
  ),
  {
    set.seed(3245)
    suppressMessages(suppressWarnings(
      gMAP(
        cbind(r, n - r) ~ 1 | study,
        data = AS,
        family = binomial,
        tau.dist = "HalfNormal",
        tau.prior = 0.5,
        beta.prior = 2,
        warmup = 500,
        iter = 1000,
        chains = 2,
        thin = 1
      )
    ))
  }
)
