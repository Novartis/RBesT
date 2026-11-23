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
`pos*` integrates θ over a design prior (probability of success, assurance).
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

`pos*` takes the same design arguments and returns a function of the **design
prior(s)**: the distribution(s) of θ the design is averaged over. A design
prior is not updated with the data of the trial part being evaluated. At an
interim it is the posterior given the data observed so far; at the planning
stage see
[borrowing-design](./borrowing-design.md).

```r
post_ctrl <- postmix(rob_ctrl, r = 8, n = 20)
post_trt <- postmix(vague_trt, r = 14, n = 20)

pos <- pos2S(vague_trt, rob_ctrl, n1 = 60, n2 = 30, decision = delta_crit)
pos(post_trt, post_ctrl)
```

Note the split: the priors passed to `pos2S` are the **analysis** priors the
final analysis will use, while the arguments to the returned function are the
**design** priors. Mixing the two up is a classic error.

Conditional power is the same call with point masses instead of posteriors —
`oc*` evaluated at the assumed θ — while predictive power integrates the
posterior, as above.

## Assurance over a design prior that links the arms

`pos2S` averages over **independent** design priors, one per arm. It cannot
express a design scenario in which the arms are linked: "control from the MAP
prior, active equal to control plus the planned effect", a design prior on a
treatment contrast (δ ~ Normal(m, s²)), or an exact null in which both arms
share the same draw (two independent draws from the same prior are not a null
effect). Use `pos2S` whenever the per-arm form is appropriate; otherwise this
is the **one permitted exception** to "RBesT computes everything": average
`oc2S` evaluations over the joint design prior by quasi-Monte Carlo (QMC).
Which scenarios to evaluate, and why, is in
[borrowing-design](./borrowing-design.md).

What is and is not allowed:

- `oc2S` supplies every probability. Its returned function is already
  integrated over the trial data, hence smooth in θ — which is why QMC
  converges fast. Do not simulate trial data and average 0/1 decisions, and do
  not hand-code the boundary or the power.
- Own code does the outer average only: the points, the inverse-CDF transform
  per margin (`qmix` for a mixture, `qnorm` etc. otherwise), the link between
  the arms, `mean()`.
- **The Monte Carlo standard error (MCSE) is mandatory.** One point set gives
  no error estimate, so the calculation must use several independently
  randomized replicates of it (16 below) and compute
  MCSE = SD of the replicate estimates / √(number of replicates). Never drop
  the replicates to save time, and never quote an estimate without its MCSE.
  If the MCSE is too large for the precision needed, raise `n_qmc` (in powers
  of two) or `n_rep` and rerun.
- The analysis priors inside `oc2S` stay untouched; the design priors only
  decide which truths are averaged. Do not use `sigma1`/`sigma2` of `oc2S` to
  inject design-prior uncertainty: they change the analysis model.
- The MCSE measures the integration error only. Do not quote more digits than
  RBesT's own numerical accuracy supports (about 1e-5), however small the MCSE.
- The design prior must have a quantile function per margin. A joint design
  prior that exists only as posterior draws is outside this recipe — escalate.

**Where the points come from.** Make an effort to use randomized QMC. Check
which of the packages `spacefillr`, `qrng` and `randtoolbox` is installed and
use the first one found. Each function below returns one replicate per `seed`
as an `n × d` matrix:

```r
## Owen-scrambled Sobol' points
qmc_spacefillr <- function(n, d, seed) {
  spacefillr::generate_sobol_owen_set(n, d, seed = seed)
}
## digitally shifted Sobol' points
qmc_qrng <- function(n, d, seed) {
  matrix(qrng::sobol(n, d, randomize = "digital.shift", seed = seed), ncol = d)
}
## Sobol' points with a random shift by hand: the scrambling of randtoolbox is
## disabled, and its sobol() leaves out the origin, which is put back here
qmc_randtoolbox <- function(n, d, seed) {
  set.seed(seed)
  points <- rbind(0, matrix(randtoolbox::sobol(n - 1, dim = d), ncol = d))
  sweep(points, 2, runif(d), "+") %% 1
}
## fallback only: plain Monte Carlo
mc_plain <- function(n, d, seed) {
  set.seed(seed)
  matrix(runif(n * d), ncol = d)
}

for (points in list(qmc_spacefillr, qmc_qrng, qmc_randtoolbox, mc_plain)) {
  stopifnot(identical(dim(points(8, 1, seed = 1)), c(8L, 1L)))
}
```

All but `qmc_spacefillr` reset the random seed of the session. If none of the
three packages is installed, ask the user to install `spacefillr` (or `qrng`).
In script-only mode ask which of them the user's platform has, and let the
script stop with an install hint if it is missing. Only if that is declined or
not possible, fall back to `mc_plain` with `n_qmc` raised to at least `2^13`,
and say in the answer that plain Monte Carlo was used. Replicates and MCSE stay
as they are. Never write your own low-discrepancy generator.

Settle with the user first: the control design prior, the planned effect and
its scale, and whether uncertainty on the effect is wanted at all.

```r
## control truth: the MAP prior itself (never the robust MAP prior)
ctrl_design <- map_ctrl
## planned effect as a log odds ratio: 25% -> 45%
delta_plan <- logit(0.45) - logit(0.25)
## optional deviation: uncertain effect with Pr(delta < 0) = 5%
delta_sd <- delta_plan / qnorm(0.95)

n_qmc <- 2^10
n_rep <- 16
qmc_points <- function(seed) qmc_spacefillr(n_qmc, 2, seed) # or the one installed

assurance_replicate <- function(seed) {
  u <- qmc_points(seed)
  theta_ctrl <- qmix(ctrl_design, u[, 1])
  theta_plan <- inv_logit(logit(theta_ctrl) + delta_plan)
  theta_unc <- inv_logit(logit(theta_ctrl) + qnorm(u[, 2], delta_plan, delta_sd))
  c(
    no_effect = mean(oc(theta_ctrl, theta_ctrl)), # both arms share one draw
    planned_effect = mean(oc(theta_plan, theta_ctrl)),
    uncertain_effect = mean(oc(theta_unc, theta_ctrl))
  )
}

replicates <- sapply(346879 + seq_len(n_rep), assurance_replicate)
stopifnot(!anyNA(replicates))

## estimate and its MCSE belong together: never report one without the other
rbind(
  estimate = rowMeans(replicates), # expected type I error, expected power
  mcse = apply(replicates, 1, sd) / sqrt(n_rep)
)
```

`oc` is the two-sample design from above; run the same average for each design
that is being compared. `planned_effect` is the Bayesian-consistent scenario
and the one to report first; `uncertain_effect` is a deviation from it and can
be left out. Other deviations change one line: another `ctrl_design`, or a
shift of the control truth, `inv_logit(logit(theta_ctrl) + shift)`. For a
normal endpoint the link is simply `theta_ctrl + delta_plan`; for a prior on
the link scale built with `family` (see above) likewise. Present the result as
a QMC (or, in the fallback, Monte Carlo) estimate based on `oc2S`, not as
`pos2S` output.

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
5. For a QMC or Monte Carlo assurance: the MCSE has been calculated and is reported, and it
   is small against the precision quoted. The no-effect result is the type I
   error *averaged over the control design prior*; with an informative
   analysis prior it need not equal the nominal level. Compare it with the
   design's OC function at a few fixed null values rather than with α.

## Next steps

[borrowing-design](./borrowing-design.md) for which designs and scenarios to
evaluate, [analysis](./analysis.md) for the posterior that feeds a PoS call,
[prior-tuning](./prior-tuning.md) for the ESS of the prior.
