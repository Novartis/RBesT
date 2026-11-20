# Design Evaluation: Decisions, Boundaries, OC and PoS

Use for success criteria, critical values, type I error and power, probability of
success, interim/futility rules, and co-data designs.

## The three layers

1. **`decision1S` / `decision2S`** build the success rule — a function of the
   posterior. Nothing is evaluated yet.
2. **`decision1S_boundary` / `decision2S_boundary`** turn that rule into the
   critical value of the observed outcome.
3. **`oc1S` / `oc2S` / `pos1S` / `pos2S`** integrate the rule over outcomes. All
   four **return a function** — you must call the result.

`oc*` conditions on an assumed true θ (frequentist power / type I error).
`pos*` integrates θ over a "knowledge" distribution (probability of success).
Same design, different question. Say which one the user actually wants.

## Building the decision

`decision1S(pc = 0.975, qc = 0, lower.tail = TRUE)` encodes
`Pr(θ ≤ qc) > pc`. With `lower.tail = FALSE` it is `Pr(θ > qc) > pc`. Nail the
direction down before writing anything.

`decision2S(pc, qc, lower.tail, link = c("identity", "logit", "log"))` does the
same for θ₁ − θ₂ (on the chosen link scale). Supplying vectors for `pc`/`qc`
creates a **composite** criterion that must hold in all components:

```r
library(RBesT)
success <- decision2S(c(0.95, 0.5), c(0, -0.10), lower.tail = FALSE)
```

`lower.tail` may have either one element or as many as `pc`. With several
elements — and especially with *mixed* elements — the rule becomes two-sided,
capturing intermediate scenarios. This is powerful and easy to get wrong; write
the resulting statement out in prose and have the user confirm it.

## One-sample design

```r
map_ctrl <- mixbeta(c(0.7, 18, 133), c(0.3, 3.5, 15))
rob_ctrl <- robustify(map_ctrl, weight = 0.2, mean = 0.5)

crit <- decision1S(0.975, 0.25, lower.tail = FALSE) # Pr(theta > 0.25) > 0.975
design <- oc1S(rob_ctrl, n = 40, decision = crit)

design(c(0.25, 0.35, 0.45)) # type I error at 0.25, power beyond
decision1S_boundary(rob_ctrl, n = 40, decision = crit) # critical number of responders
```

## Two-sample design

```r
vague_trt <- mixbeta(c(1, 1, 1))
delta_crit <- decision2S(0.95, 0, lower.tail = FALSE)

oc <- oc2S(vague_trt, rob_ctrl, n1 = 60, n2 = 30, decision = delta_crit)
oc(theta1 = c(0.4, 0.5), theta2 = c(0.25, 0.25))
```

`decision2S_boundary(prior1, prior2, n1, n2, decision)` gives the critical value
of arm 1 as a function of the arm 2 outcome.

## Normal endpoints and other exponential families

Normal methods take `sigma` (`sigma1`/`sigma2` in the two-sample case). It
defaults to the prior's reference scale, so pass it explicitly whenever the
prior carries none.

Since 1.10-0 these functions also accept a `family`, which lets the normal
approximation stand in for any one-parameter exponential family — together with
`offset1`/`offset2` on the link scale. The recurrent-event design is the
canonical case:

```r
lambda_ctrl <- 1.8
kappa_true <- 1.9
followup <- 1 * 0.85 + 0.5 * 0.15
log_exposure <- log(followup)
sigma_ctrl <- round(
  sqrt((1 + kappa_true * lambda_ctrl * followup) / (lambda_ctrl * followup)), 2
)

uninf_ctrl <- mixnorm(c(1, log(lambda_ctrl), 10), sigma = sigma_ctrl)
uninf_trt <- mixnorm(c(1, log(lambda_ctrl), 10), sigma = sigma_ctrl)
nb_family <- MASS::negative.binomial(theta = 1 / kappa_true)

oc_nb <- oc2S(uninf_trt, uninf_ctrl,
  n1 = 225, n2 = 225,
  decision = decision2S(0.975, 0, lower.tail = TRUE),
  family = nb_family,
  offset1 = log_exposure, offset2 = log_exposure
)
oc_nb(log(lambda_ctrl) + c(0, -0.3), rep(log(lambda_ctrl), 2))
```

For time-to-event analysed by the normal approximation of the log-rank statistic,
the standard error of the log hazard ratio scales as `2/sqrt(N_events)` under 1:1
randomization — hence a sampling standard deviation of **σ = 2**, with `n` counted
in events.

## Probability of success

`pos*` takes the same design arguments and returns a function of the
distribution(s) describing current knowledge about θ. At an interim, that
knowledge distribution is the posterior given the data observed so far.

```r
post_ctrl <- postmix(rob_ctrl, r = 8, n = 20)
post_trt <- postmix(vague_trt, r = 14, n = 20)

pos <- pos2S(vague_trt, rob_ctrl, n1 = 60, n2 = 30, decision = delta_crit)
pos(post_trt, post_ctrl)
```

Note the split: the priors passed to `pos2S` are the **analysis** priors the
final analysis will use, while the arguments to the returned function are the
current **knowledge** distributions. Mixing the two up is a classic error.

Conditional power is the same call with point masses instead of posteriors —
`oc*` evaluated at the assumed θ — while predictive power integrates the
posterior, as above.

## Co-data and the MAC approach

Two equivalent routes when concurrent trials contribute data:

- **MAP (sequential)**: derive the prior from the historical data, then
  `postmix` it with the current trial.
- **MAC (combined)**: refit everything at once with
  `update(map_fit, data = all_data)`, which re-runs the MCMC.

They agree. MAC is the natural choice when co-data arrive as whole extra
studies; MAP when they arrive as summary statistics for the same group.

## Cost and tuning

Boundary construction has been adaptive since 1.11-0 and is no longer the
bottleneck it once was — roughly 10× faster for normal endpoints and far more
for near-flat discrete boundaries, with `oc2S`/`pos2S` 2–3× faster. Large
scenario grids are therefore usually affordable; measure before optimizing.

- `Ngrid` sets the finest resolution of the adaptive scheme; the default dropped
  from 10 to 5 in 1.11-0 for normal endpoints. Pass `Ngrid = 10` for the old
  density. Differences are of order 1e-5.
- `options(RBesT.decision2S_boundary = "grid")` restores the legacy uniform
  sweep, for reproducing pre-1.11 numbers exactly.

## Sanity checks

1. Evaluate the OC function at the null — does it return the intended α?
2. Check the boundary value against intuition; an implausible critical value
   usually means `lower.tail` is inverted.
3. For normal endpoints, confirm σ and `n` units (patients vs events).
4. State plainly whether a number is conditional power, predictive power, or
   unconditional PoS.

## Next steps

[analysis](./analysis.md) for the posterior that feeds a PoS call,
[prior-tuning](./prior-tuning.md) for the ESS the design borrows.
