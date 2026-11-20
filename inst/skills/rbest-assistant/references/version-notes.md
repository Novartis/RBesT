# Version Notes

The rest of this skill is written against **RBesT 1.12-1**. Have the user run
`packageVersion("RBesT")` before giving version-dependent advice.

## Why did my numbers change?

| Symptom | Cause | Since |
|---|---|---|
| `Rhat`, `n_eff` and posterior summaries differ from an older run | summaries now come from the modern `posterior` functions, including a redefined `Rhat` | 1.10-0 |
| `oc2S`/`pos2S` results differ in the 5th decimal | adaptive decision boundaries and `Ngrid` default 10 → 5 | 1.11-0 |
| `oc2S`/`pmixdiff`/`ess` differ slightly | Gaussian quadrature replaced adaptive Gauss-Kronrod on identity-link paths | 1.10-0 |
| A `gMAP` fit explores differently, or divergences appear/disappear | sum-to-zero parametrization, automatic partial centering, `adapt_delta` 0.99 → 0.95 | 1.12-0 |
| `decision2S_boundary`/`oc2S`/`pos2S` results for a normal endpoint differ, possibly a lot (e.g. an OC of `0` that is now `1`) | two bugs fixed for normal endpoints: outcomes where the decision does not depend on sample 1 were dropped (with `family`), and the adaptive boundary could be traced too coarsely for mixture priors | 1.12-1 |
| Normal / MVN mixture fits are numerically different from pre-1.7-2 | the EM is now deterministic and no longer seed-dependent | 1.7-2 |

## Why did my code break?

| Error | Cause | Fix |
|---|---|---|
| `$fit` is `NULL`, or `rstan::extract` on a `gMAP` object fails | the raw `stanfit` slot was removed when posterior storage moved to `posterior` draws objects | use `as_draws()`, `as_draws_matrix()`, … or `as.matrix()` |
| Warmup draws are missing | warmup is no longer stored | `options(RBesT.MC.save_warmup = TRUE)`, then `inc_warmup = TRUE` |
| Code indexing draws by column position picks the wrong variable | variable order is now `theta`, `tau`, `beta`, `theta_pred`, `theta_resp_pred`, `lp__` | index by name |
| `write_mix_json` now errors or warns | it refuses to silently round all weights to zero | raise `digits` |

## Restoring previous behaviour

```r
options(
  RBesT.MC.s2z = FALSE, # pre-1.12 parametrization (also restores adapt_delta 0.99)
  RBesT.MC.ncp = 1, # pre-1.12 non-centered parametrization
  RBesT.MC.save_warmup = TRUE, # pre-1.9 warmup storage
  RBesT.integrate_method = "adaptive", # pre-1.10 integration
  RBesT.decision2S_boundary = "grid" # pre-1.11 boundary sweep
)
```

Plus `Ngrid = 10` as an argument to `oc2S`/`pos2S`/`decision2S_boundary` for the
pre-1.11 grid density. Use these to reproduce old numbers, not as defaults.

## What became available when

**1.12-1** — bugfix release for two-sample designs with normal endpoints.
`decision2S_boundary`, `oc2S` and `pos2S` were wrong, possibly entirely, when
used with `family` and the decision does not depend on the outcome of sample 1
for some outcomes of sample 2; and with mixture priors under the default
adaptive boundary, results could depend on earlier calls to the design object.
`oc2S`/`pos2S` now warn when the boundary changes a second time at extreme
outcomes with non-negligible probability. Anyone who computed two-sample
normal-endpoint designs with 1.12-0 (or earlier with `family`) should rerun
them on 1.12-1.

**1.12-0** — sum-to-zero parametrization of group effects (intercept + normal
random effects + single τ stratum); automatic partial centering
(`RBesT.MC.ncp = 3`) with a warned fallback; quadrature-based scaling of
`log(tau)`, fixed and group effects; `adapt_delta` default 0.95; reordered draws
variables; `write_mix_json` precision guard; webR/WebAssembly build.

**1.11-0** — adaptive decision boundaries (≈10× faster for normal endpoints,
orders of magnitude for near-flat discrete ones; `oc2S`/`pos2S` 2–3× faster);
shape-preserving monotone interpolation and exact linear boundaries for
single-component priors with constant σ; `Ngrid` default 5.

**1.10-0** — `family` argument on `oc1S`, `oc2S`, `pos1S`, `pos2S`,
`decision1S_boundary`, `decision2S_boundary` for the normal case, enabling
arbitrary exponential-family endpoints under the normal approximation; the
negative-binomial futility-interim article; Gaussian quadrature integration;
modern `posterior` summaries; `gMAP` posterior storage refactor (`fit` slot
removed).

**1.9-0** — vector `lower.tail` with as many elements as `pc`, including mixed
elements, for both `decision1S` and `decision2S` (two-sided decisions);
`RBesT.MC.save_warmup = FALSE` default; experimental `as_draws*()`.

**1.8-2** — experimental `write_mix_json` / `read_mix_json`; `"msr"`
parametrization for `mixmvnorm` with the new `msr2mvnorm`; zero standard
deviations allowed in multivariate normal components.

**1.8-0** — `family` argument for `ess()` on normal mixtures, giving the ESS of a
mixture that represents a transformed (e.g. logit) scale.

**1.7-4** — `mixstanvar` generates distribution functions, so truncated mixture
priors work in brms.

**1.7-1** — `mixstanvar` introduced; named dimensions for multivariate normal
mixtures; experimental MVN EM diagnostic plots.

## Older ESS fixes

ELIR computation was stabilized in 1.7-4 (zero-weight components, boundary
aborts) and again in 1.8-0 (beta mixture information function). An ESS number
from before those releases may not reproduce.
