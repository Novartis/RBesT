# Limitations and Escalation

## State these proactively

1. `gMAP` runs Stan MCMC. It is not instant, and `update()` re-runs it in full.
2. Reproducibility requires an explicit `set.seed()` plus fixed iteration and
   chain settings.
3. Convergence must be checked before anything is interpreted: divergences,
   `Rhat`, and the τ posterior against its prior.
4. ESS is method-dependent. `elir` is the default and the recommendation;
   `moment` is conservative; `morita` can be optimistic with many components.
   Never quote an ESS without its method.
5. The robustification weight is a modelling choice with real consequences, not a
   tuning parameter to pick quietly.
6. `automixfit` selects components by AIC and can still fit poorly — the density
   overlay in `plot(...)$mix` is the check, and `mixfit(..., Nc = )` the override.
7. Normal `gMAP` without `weights` yields a mixture with no reference σ, so every
   downstream call needs `sigma` supplied explicitly.
8. A vector `lower.tail` in `decision1S`/`decision2S` creates a two-sided rule.
   It is easy to write and easy to misread — restate it in prose.
9. The variance MAP workflow rests on a normal approximation to a log-Gamma
   variate: useful from about 10 degrees of freedom, good beyond 20.
10. Operating-characteristic grids are much cheaper since 1.11-0, but a large
    two-sample grid over mixture priors is still real computation. Measure a
    small grid first.
11. For a confirmatory trial that borrows historical data, the prior, the
    robust weight and the evaluation of the design have to be agreed with the
    health authority. The skill helps to prepare that discussion; it does not
    settle it.

## Escalate to the package authors

Contact the maintainer, Sebastian Weber (<sebastian.weber@novartis.com>), by
e-mail when any of the following holds:

1. The endpoint, likelihood or link function is unsupported or unclear.
2. Borrowing assumptions are ambiguous and consequential.
3. The hierarchy or covariate structure goes beyond the documented workflows.
4. The decision framework does not map onto `decision1S`/`decision2S`.
5. Runtime or numerical problems block practical execution.
6. Prior-data conflict strategy is still unresolved after the follow-ups.
7. The prior has to be embedded in an external Bayesian model beyond the
   documented `mixstanvar`/brms pattern.
8. Multivariate or multi-endpoint borrowing is required — `mixmvnorm` exists but
   the surrounding workflows are thin.
9. Time-to-event work needs a survival model beyond the log-rank normal
   approximation.
10. Sample-size re-estimation or adaptive designs beyond PoS-based futility.
11. The question is about simulation-based calibration (`inst/sbc/`), which is
    development infrastructure rather than a user-facing feature.
12. The `tau.strata` assignment for a specific constellation of studies is not
    obvious.
13. A result seems to need numerical code beyond RBesT's own functions and the
    permitted exceptions listed in the skill — including an assurance where
    the joint design prior exists only as posterior draws. Do not write that
    code; ask.
14. The historical data show a convincing time trend, or the user wants to
    adjust for baseline characteristics beyond a study-level covariate in the
    `gMAP` formula. These are advanced topics: point to the case studies at
    <https://opensource.nibr.com/bamdd/> first, where such models are demonstrated.

When escalating, ask the user to prepare a minimal reproducible example with the
data schema, the intended decision rule, and their current RBesT code.
