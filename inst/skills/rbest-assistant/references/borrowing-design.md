# Planning a Trial that Borrows Historical Control Data

Use when a two-arm trial is being planned with an informative (MAP) prior for
the control arm: what to aim for, which designs to compare, and under which
scenarios to evaluate them. The mechanics of `decision2S`, `oc2S` and `pos2S`
are in [design](./design.md); deriving and robustifying the prior in
[map-prior](./map-prior.md) and [prior-tuning](./prior-tuning.md).

This file gives a direction of thinking, not a fixed procedure. The user and
you are expected to explore the design space along these lines.

## Follow-ups to settle first

1. **The goal.** Either *reduce the total sample size* by enrolling fewer
   controls, or *keep the total sample size N fixed* and shift the allocation
   towards active (2:1, say), with the prior supplementing the smaller control
   arm. The second can give more information — and gives more data on the
   active treatment — from the same resources.
2. **The setting.** Exploratory (the sponsor's own risk) or confirmatory. For a
   confirmatory trial the prior, the robust weight and the evaluation are to be
   agreed with the health authority early; say so.
3. **The planned treatment effect and its scale** (risk difference, log odds
   ratio, mean difference, log hazard ratio), as used for the sample size.
4. **The candidate allocation(s)** and total sample size(s).

## Analysis priors and design priors

- An **analysis prior** enters the analysis of the trial: it is updated with
  the trial data. Here: the robust MAP prior for control, a weakly informative
  prior for active.
- A **design prior** describes which truths the design is evaluated under. It
  is never updated with data and does not change the decision rule; it only
  weights the operating characteristics.

The FDA draft guidance "Use of Bayesian Methodology in Clinical Trials of Drug
and Biological Products" (January 2026,
<https://www.fda.gov/regulatory-information/search-fda-guidance-documents/use-bayesian-methodology-clinical-trials-drug-and-biological-products>)
makes this distinction. Among other approaches it distinguishes success
criteria *calibrated* to the type I error rate from criteria interpreted
directly as posterior probabilities and evaluated under design priors.
Borrowing through a MAP prior is naturally read in the second way. Type I error
and power at fixed parameter values are still reported, but an informative
prior cannot control the type I error at every control value; whether that is
acceptable is to be agreed with the agency. Do not claim otherwise, do not
present the guidance as prescribing one approach, and say that it is a draft
and may change.

## Designs to compare

Evaluate the planned design next to these references. All use the same decision
rule; only the control analysis prior and the allocation change.

Evaluate with the prior **as it is pre-specified**: the mixture read back with
`read_mix_json()` from the file that documents it, rounding included — not the
`automixfit` output it was written from and not the MCMC sample. That restored
object is what the final analysis will use.

| Control analysis prior | Allocation | What it shows |
|---|:---:|---|
| Non-informative | 1:1 | The conventional design |
| Non-informative | R:1 | The trial data alone — what is left if nothing is borrowed |
| MAP | R:1 | Full use of the historical information |
| Robust MAP | R:1 | The planned design |

For the fixed-N goal all four have the same total N; for the reduction goal the
1:1 design has the conventional, larger N.

```r
library(RBesT)
map_ctrl <- mixbeta(c(0.7, 18, 133), c(0.3, 3.5, 15)) # stands in for a MAP prior
rob_ctrl <- robustify(map_ctrl, weight = 0.2, mean = 0.5)
flat <- mixbeta(c(1, 1, 1))
success <- decision2S(0.975, 0, lower.tail = FALSE)

designs <- list( # fixed N = 90
  flat_1to1 = oc2S(flat, flat, 45, 45, success),
  flat_2to1 = oc2S(flat, flat, 60, 30, success),
  map_2to1 = oc2S(flat, map_ctrl, 60, 30, success),
  rmap_2to1 = oc2S(flat, rob_ctrl, 60, 30, success)
)

theta_ctrl <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.40) # around and away from the prior
theta_trt <- inv_logit(logit(theta_ctrl) + log(3)) # planned effect: odds ratio 3

## each design is a function of the true (active, control) parameters
data.frame(theta_ctrl, sapply(designs, do.call, list(theta_ctrl, theta_ctrl))) # type I error
data.frame(theta_ctrl, sapply(designs, do.call, list(theta_trt, theta_ctrl))) # power
```

Plot or tabulate over the whole plausible range of the control parameter, not
at the prior mean only. Take the range from what the historical data support —
the central 95% interval of the MAP prior (`summary(map_ctrl)`), somewhat
widened — and say so; far outside it the type I error of a borrowing design
looks alarming, but those control values contradict everything observed so far.

How to read the curves: an informative prior pulls the control estimate towards
the historical mean. Where the true control value lies on the side on which
this makes the treatment look better, the type I error rises above the nominal
level; on the other side it falls below, and power moves the same way. This is
the price of borrowing, not an error in the code: a gain in power from an
informative prior is not possible under strict type I error control.

## What the prior contributes depends on the data

The ESS (`ess()`, ELIR) describes how informative the prior is, and is to be
reported. It is **not** a number of patients that can be taken out of, or added
to, a frequentist power calculation. Never size a trial as "N minus ESS" and
never describe a design as worth "N plus ESS" patients.

The reason: how much the prior adds to the posterior depends on the control
data the trial will observe. If they agree with the prior, the gain is largest.
The more they disagree, the smaller it gets, down to a posterior that is less
clear than one without the prior, and a biased one: robustification limits the
bias but does not remove it. The table above shows exactly this as a function
of the true control value, and the R:1 design without borrowing shows the price
of the unequal allocation when the prior does not help.

To make this visible, compare the ESS of the control posterior with that of the
prior for every possible control outcome, next to how likely each outcome is if
the MAP prior holds:

```r
n_ctrl <- 30
r_ctrl <- 0:n_ctrl
ess_gain <- sapply(r_ctrl, function(r) ess(postmix(rob_ctrl, r = r, n = n_ctrl))) - ess(rob_ctrl)
likely <- dmix(preddist(map_ctrl, n = n_ctrl), r_ctrl) # prior predictive of the outcome
round(data.frame(r_ctrl, ess_gain, likely), 2)
```

Outcomes in line with the prior add more information than the `n_ctrl` patients
enrolled; outcomes in conflict add less, and can leave less than the prior
alone had.

The information in the prior is fixed, while that in the data grows with the
sample size. An analysis done earlier than planned — an interim, or a trial
stopped short — therefore gives the prior more influence, a later one less.
Evaluate interim analyses with the sample size they will actually have.

## Scenarios to evaluate

Fixed parameter values, as above, answer "what if the truth is here". Which
truths to expect is a separate question:

1. **Bayesian-consistent scenario — start here.** If the MAP prior is believed,
   it is also the prediction of the control parameter in the new trial. So let
   the control truth follow the MAP prior and set active to control plus the
   planned effect. That effect is a point value: a guess taken from the sample
   size rationale, with no reliable uncertainty attached — say so rather than
   dressing it up. The no-effect case uses the same control draw for both arms.
   The averages over this scenario are the expected type I error and the
   expected power (assurance) of each design; their differences between designs
   are the expected loss and gain from borrowing and from the allocation.
   Expected power is usually lower than the power at the planned alternative,
   because the uncertainty about the truth is averaged in. That is the intended
   result, not a sign that something went wrong.
2. **Deviations from it** are what other design priors are for. Offer them as
   sensitivity scenarios, one at a time:
   - a MAP prior derived with another heterogeneity prior
     (`update(map_fit, tau.prior = ...)`, then `automixfit`);
   - uncertainty on the treatment effect instead of a point value;
   - a shift of the control truth, e.g. motivated by a time trend (below).

The **robust MAP prior is not a design prior**: robustification acts through
updating with data, which never happens to a design prior. Use the MAP prior.

`pos2S` cannot compute these averages, because it treats the two arms as
independent. Use the recipe "Assurance over a design prior that links the arms"
in [design](./design.md); it covers all scenarios listed here.

## Time trend in the historical data

Standard of care changes. Ask for the year of each historical study, sort the
data by it before fitting (the forest plot follows the order of the data), and
look at `forest_plot(map_fit)` or `plot(map_fit)$forest_model`. If the control
outcome drifts in one direction, define a design scenario in which the control
truth is moved further that way: fixed control values beyond the latest studies
in the table above, or a shift on the link scale in the assurance recipe.

With few studies a trend is weak evidence: treat it as a "what if" scenario,
not as an estimate. A convincing trend questions the exchangeability
assumption behind the analysis prior itself. What can be done about it:

- Discounting older studies more strongly with `tau.strata`, see
  [map-prior](./map-prior.md).
- A study-level covariate in the `gMAP` formula (meta-regression), see
  [map-prior](./map-prior.md) for what this covers and its limits.
- Beyond that — a trend or random-walk model, adjustment for patient baseline
  characteristics, subgroup results reported by only some studies — these are
  advanced topics. Say that they exist, recommend the worked case studies at
  <https://opensource.nibr.com/bamdd/>, and [escalate](./escalation.md) instead of building such a model
  yourself.

## Sanity checks

1. Every design in the comparison uses the same decision rule and, for the
   fixed-N goal, the same total N.
2. Type I error and power are shown across the plausible control range, and the
   range in which the informative design exceeds the nominal level is stated.
3. Every averaged result names its design prior and carries its Monte Carlo
   standard error.
4. No sentence equates the ESS with patients saved or gained.

## Next steps

[design](./design.md) for the assurance recipe, boundaries and interim PoS;
[prior-tuning](./prior-tuning.md) for the robust weight and the ESS.
