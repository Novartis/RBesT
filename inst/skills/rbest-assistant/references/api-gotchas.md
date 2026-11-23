# API Gotchas and Global Options

Load this instead of reconstructing a call from memory. Everything here was
checked against the RBesT 1.12-1 sources.

## Verified signatures

```text
gMAP(formula, family = gaussian, data, weights, offset, tau.strata,
     tau.dist = c("HalfNormal", "TruncNormal", "Uniform", "Gamma", "InvGamma",
                  "LogNormal", "TruncCauchy", "Exp", "Fixed"),
     tau.prior, tau.strata.pred = 1, beta.prior, prior_PD = FALSE,
     REdist = c("normal", "t"), t.df = 5, contrasts = NULL,
     iter, warmup, thin, init, chains, cores)      # last six default from options

mixbeta(..., param = c("ab", "ms", "mn"))
mixnorm(..., sigma, param = c("ms", "mn"))
mixgamma(..., param = c("ab", "ms", "mn"), likelihood = c("poisson", "exp"))
mixmvnorm(..., sigma, param = c("ms", "mn", "msr"))

mixfit(sample, type, Nc, ...)
automixfit(sample, Nc = seq(1, 4), k = 6, thresh = -Inf, verbose = FALSE, ...)
mixcombine(..., weight, rescale = TRUE)
robustify(priormix, weight, mean, n = 1, ...)          # normMix also takes sigma
ess(mix, method = c("elir", "moment", "morita"), ...)  # normMix also takes family, sigma

postmix(priormix, data, ...)                           # betaMix: n, r
                                                       # normMix: n, m, se
                                                       # gammaMix: n, m
preddist(mix, n = 1, ...)                              # normMix also takes sigma
pmixdiff(mix1, mix2, q, lower.tail = TRUE)             # also d/q/rmixdiff

decision1S(pc = 0.975, qc = 0, lower.tail = TRUE)
decision2S(pc = 0.975, qc = 0, lower.tail = TRUE,
           link = c("identity", "logit", "log"))
decision1S_boundary(prior, n, decision, ...)
decision2S_boundary(prior1, prior2, n1, n2, decision, ...)
oc1S(prior, n, decision, ...)                          # normMix: sigma, family, offset
oc2S(prior1, prior2, n1, n2, decision, ...)            # normMix: sigma1, sigma2, eps,
                                                       #   Ngrid = 5, family, offset1, offset2
pos1S(prior, n, decision, ...)                         # same extras as oc1S
pos2S(prior1, prior2, n1, n2, decision, ...)           # same extras as oc2S

mixstanvar(..., verbose = FALSE)
write_mix_json(mix, con, ...)
read_mix_json(con, ..., rescale = TRUE)
```

## Things that are commonly gotten wrong

| Trap | Reality |
|---|---|
| `plot(gmap_fit)$tau_forest` | No such element. The list is `traceBeta`, `traceTau`, `traceLogTau`, `densityBeta`, `densityTau`, `densityLogTau`, `densityThetaStar`, `densityThetaStarLink`, `forest_model` |
| `oc1S(...)` returns numbers | It returns a **function**; call it with θ |
| `pos2S(...)` arguments | Priors passed to `pos2S` are the analysis priors; the returned function takes the design priors (at an interim: the current posteriors), independent per arm |
| `sigma` must be global | It is a formal argument of the normal `oc*`/`pos*`/`preddist`/`decision*_boundary` methods, defaulting to `sigma(prior)` |
| `update()` is cheap | `stats::update.default` re-runs the full MCMC from the stored `$call`; there is no `update.gMAP` method |
| `tau.prior = c(a, b)` means two strata | With a *single* stratum it is read as the parameter pair (a, b). Multiple strata need one entry per stratum |
| normal `gMAP` without `weights` | Allowed, but the fit cannot estimate a reference σ, so the resulting mixture has none and every downstream call must be given `sigma` |
| `mixfit` needs raw samples | It also accepts `gMAP` and prediction objects directly |
| `ess()` number is comparable across methods | It is not — `elir` (default), `moment` and `morita` answer slightly different questions. Always state the method |
| `lower.tail` is cosmetic | It flips the inequality. Vector `lower.tail` with mixed elements makes the decision two-sided |

Exported but superseded: `oc1Sdecision()` and `oc2Sdecision()`. If they appear in
a user's legacy code, translate to `decision1S()` / `decision2S()`.

`logit()`, `inv_logit()`, `lodds()` and `BinaryExactCI()` are exported and often
what someone is hand-rolling.

## Self-check

```r
library(RBesT)
stopifnot(
  is.function(oc1S(mixbeta(c(1, 1, 1)), 20, decision1S(0.95, 0.3))),
  formals(robustify)$n == 1,
  identical(eval(formals(ess)$method), c("elir", "moment", "morita"))
)
```

## Global options

`?RBesT-package` is the authoritative list. The ones that come up:

| Option | Default | Effect |
|---|---|---|
| `RBesT.MC.warmup` | 2000 | MCMC warmup iterations |
| `RBesT.MC.iter` | 6000 | total MCMC iterations |
| `RBesT.MC.chains` | 4 | chains |
| `RBesT.MC.thin` | 4 | thinning |
| `RBesT.MC.save_warmup` | `FALSE` | keep warmup draws in the fit |
| `RBesT.MC.control` | `adapt_delta = 0.95` | passed to Stan; raise to 0.99 on divergences |
| `RBesT.MC.s2z` | `TRUE` | sum-to-zero parametrization of group effects |
| `RBesT.MC.ncp` | 3 | group-effect parametrization: 0 centered, 1 non-centered, 2/3 partial |
| `RBesT.MC.init` | 1 | range of the uniform initial values |
| `RBesT.MC.rescale` | `TRUE` | automatic rescaling of raw parameters |
| `RBesT.integrate_method` | quadrature | `"adaptive"` restores Gauss-Kronrod integration |
| `RBesT.GQ_nodes` | — | starting node count for Gaussian quadrature |
| `RBesT.GQ_abs_tol`, `RBesT.GQ_rel_tol` | — | quadrature refinement tolerances |
| `RBesT.decision2S_boundary` | adaptive | `"grid"` restores the pre-1.11 uniform sweep |

Since 1.10-0 numerical integration against mixture densities uses Gaussian
quadrature on identity-link paths — Gauss-Hermite for normal, Gauss-Jacobi for
beta, Gauss-Laguerre for gamma — in `oc2S`, `pos2S`, `dmixdiff`, `pmixdiff` and
`ess`. Node counts refine automatically until successive estimates agree within
`max(abs_tol, rel_tol * |I|)`, falling back to adaptive integration at the cap.
Non-identity links and beta-mixture ESS still use adaptive integration.

## Next steps

[version-notes](./version-notes.md) for what changed and when.
