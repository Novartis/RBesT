## Verbose S2Z fit at the normal endpoint, used by the zeta invariant test.
##
## The normal endpoint is the one where a dropped tau-dependent normalising
## constant would be least visible in the binomial fits, so zeta is checked
## here as well.
fixture <- withr::with_options(
  list(
    RBesT.verbose = TRUE,
    RBesT.MC.ncp = 1,
    RBesT.MC.save_warmup = FALSE
  ),
  {
    crohn_se <- crohn
    crohn_se$y.se <- 88 / sqrt(crohn_se$n)
    set.seed(46711)
    suppressMessages(suppressWarnings(
      gMAP(
        cbind(y, y.se) ~ 1 | study,
        data = crohn_se,
        family = gaussian,
        tau.dist = "HalfNormal",
        tau.prior = 44,
        beta.prior = 88,
        warmup = 1000,
        iter = 2000,
        chains = 4,
        thin = 1
      )
    ))
  }
)
