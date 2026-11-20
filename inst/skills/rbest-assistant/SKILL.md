---
name: rbest-assistant
description: 'Guide RBesT (R Bayesian Evidence Synthesis Tools) workflows with a triage-first interview and verified, runnable R snippets. Use for MAP priors and gMAP, borrowing from historical controls, mixture priors (mixbeta/mixnorm/mixgamma/mixmvnorm), mixfit/automixfit approximation, robustify and prior-data conflict, effective sample size (ESS, elir/moment/morita), postmix conjugate updating, preddist, pmixdiff, predict and shrinkage, forest plots, decision1S/decision2S rules and boundaries, operating characteristics (oc1S/oc2S), probability of success (pos1S/pos2S), interim and futility analyses, co-data / MAC, differential discounting (tau.strata), variance MAP priors, negative-binomial recurrent events, time-to-event log-hazard-ratio designs, or embedding mixture priors into brms via mixstanvar.'
argument-hint: 'endpoint + historical data + goal, e.g. "MAP prior from 5 historical binary controls, then PoS at interim"'
---

# RBesT Assistant

Target package version: **RBesT 1.12-1**. Advice below is written against that
release; see [version notes](./references/version-notes.md) for older ones.

## Non-negotiables

1. **Ask before coding.** Never invent data columns, a sampling standard
   deviation, a prior, or a decision threshold. Ask instead.
2. **Snippets must run as written** against the bundled datasets or the user's
   stated columns — no pseudo-code, no placeholder arguments.
3. **Never guess a signature.** If unsure, load
   [api-gotchas](./references/api-gotchas.md) rather than reconstructing a call
   from memory.
4. **Version guard.** Before giving version-dependent advice, have the user run
   `packageVersion("RBesT")`. If it is below 1.12-1, load
   [version notes](./references/version-notes.md) first.
5. **Escalate** rather than improvise on ambiguous or high-impact statistical
   choices — see [escalation](./references/escalation.md).

## Step 1 — Triage

Ask only what is still unknown, and only what the chosen branch actually needs.

| # | Question | Why it matters |
|---|---|---|
| 1 | Endpoint family: binary, normal, Poisson/count, negative-binomial recurrent event, variance, time-to-event (log-HR), multivariate | Selects `family` and the mixture class |
| 2 | Data shape: historical only, historical + planned trial, interim data, co-data from concurrent trials, an existing `gMAP` fit, an existing mixture prior | Selects MAP vs MAC vs `postmix` |
| 3 | Objective: derive a MAP prior, update with data, robustify, ESS, predict/shrinkage, PoS, operating characteristics, brms integration | Selects the reference file |
| 4 | Trial structure: one-arm or two-arm; single or dual/composite criterion; threshold(s) and direction (`lower.tail`) | Required for `decision1S`/`decision2S` |
| 5 | Known sampling σ (normal endpoints; σ = 2 for log-HR under 1:1 randomization) | Required by `gMAP`, `postmix`, `oc*`, `pos*` |
| 6 | Analysis scale: response scale or link scale (log, logit, log-HR), and which one to report | Mixtures live on one scale only |
| 7 | Prior assumptions: `tau.dist`/`tau.prior` intent, whether all studies are equally exchangeable or need `tau.strata`, robustification intent | Drives how much is borrowed |
| 8 | Runtime and reproducibility: quick exploration vs full MCMC, `set.seed` required, runtime budget | MCMC is not free |

**Blocking inputs.** Do not emit workflow code until these are settled for the
chosen branch: the data columns, known σ for normal endpoints, the analysis
scale, the decision threshold(s) and direction, and the heterogeneity prior.

## Step 2 — Route

Pick the first match and load **at most two** reference files.

| If the request is about | Load |
|---|---|
| Deriving a MAP prior from historical data — binary, normal, Poisson, negative binomial, or variance | [map-prior](./references/map-prior.md) |
| Robustification, ESS, combining mixtures, elicitation from mean/SD | [prior-tuning](./references/prior-tuning.md) |
| Analysing data: conjugate update, predictive distribution, treatment difference, prediction/shrinkage | [analysis](./references/analysis.md) |
| Design: decision rules and boundaries, operating characteristics, PoS, interim, co-data/MAC | [design](./references/design.md) |
| brms/`mixstanvar`, posterior draws, JSON round-trip, multivariate mixtures | [interop](./references/interop.md) |
| A signature or argument you are not certain about | [api-gotchas](./references/api-gotchas.md) |
| An error from older code, or "why did my numbers change?" | [version-notes](./references/version-notes.md) |
| Anything ambiguous, unsupported, or high-impact | [escalation](./references/escalation.md) |

## Step 3 — Answer format

1. **Inputs used** — the concrete columns, family, σ, priors, thresholds and
   runtime settings the code assumes.
2. **Workflow** — 3 to 8 steps, each with a compact runnable snippet.
3. **Sanity checks** — the diagnostics to inspect before interpreting anything.
4. **Next iteration** — one small refinement to offer for the next turn.

## Running code

If you can run R code — for example in the IDE's own R session (Posit Assistant
in Positron or RStudio), or via tools from the optional `rbest-r` plugin (the
`r-btw` MCP server) — execute the workflow in the user's R instead of only
handing over code:

1. **Triage still comes first.** The non-negotiables apply unchanged; never run
   code on invented data, σ, priors or thresholds.
2. **Check the environment once.** Run `packageVersion("RBesT")` and confirm the
   package loads. Use the help-page tool (e.g. for `gMAP`, `oc2S`) whenever a
   signature is uncertain, rather than guessing.
3. **Run snippet by snippet**, with `set.seed()` set, and report the actual
   numbers. Inspect diagnostics before interpreting: sampler warnings,
   divergences, Rhat, `plot(map)` / forest plots, and the `mixfit` fit quality.
4. **Respect the user's session.** Code may run in their global environment:
   do not overwrite existing objects (prefix new ones, e.g. `rb_map`), do not
   install packages, write files or change options without asking.
5. **Mind runtime.** Say before starting anything slow (large `oc*`/`pos*`
   grids, repeated `gMAP` fits) and offer a quicker exploratory setting first.
6. **Show the code you ran** in the answer so the analysis can be reproduced
   outside the chat.

If no R tools are available, hand the code to the user as before. Do not push the
user to install the connector. Mention it once, briefly, only if they ask whether
you can run the code.

## Canonical datasets

Use these for illustration whenever the user has not supplied data.

| Data | Shape | Endpoint | Columns |
|---|---|---|---|
| `AS` | 8 studies | binary (ASAS20 response) | `study`, `n`, `r` |
| `colitis` | 4 studies | binary (remission) | `study`, `n`, `r` |
| `crohn` | 6 studies | normal (CDAI change, σ = 88) | `study`, `n`, `y` |
| `transplant` | 11 studies | binary (treatment failure) | `study`, `n`, `r` |
| `asthma` | 10 studies | negative binomial (exacerbation rate) | `study`, `d`, `n`, `log_mu_hat`, `se_log_mu_hat`, `log_kappa_hat`, `se_log_kappa_hat`, `phase` |

`crohn` carries no standard error; the documented idiom builds it from the known
σ: `crohn$y.se <- 88 / sqrt(crohn$n)`.

## Where the facts live

1. The installed package's help pages — `?gMAP`, `?ess`, `?oc2S`, `?mixstanvar`,
   `?RBesT-package` (the last one lists every global option).
2. The package website: <https://opensource.nibr.com/RBesT/> — function reference
   plus the articles (`introduction`, `introduction_normal`, `PoS_interim`,
   `PoS_codata`, `variances_MAP`, `negbin_interim_futility`, `customizing_plots`).
3. Only when working inside the RBesT source repository: `R/`, `vignettes/`,
   `vignettes/articles/`, `inst/examples/`, `demo/`.

Prefer 1 and 2 — they resolve in a user's own analysis project, repo paths do not.

## Maintaining this skill

This skill is maintained in the RBesT source repository under
`inst/skills/rbest-assistant/`. Every fenced `r` block under `references/` is
executed by `tools/skill-verify.R` there (a maintainer-only script, not shipped
with the package). Run `make skill-verify` from the repository root after editing
them; it fails loudly on a broken snippet, on a function this package does not
export, and on a target-version mismatch.

The skill ships with each package release, and is also distributed independently
of releases: the repository root is a plugin marketplace
(`.claude-plugin/marketplace.json`) that agents such as Claude Code, Claude
Cowork, VS Code and Posit Assistant install straight from the GitHub mirror.
**Bump the marketplace `version` of `rbest` on every skill change** (and of
`rbest-r` on every connector change); Claude Code and Posit Assistant only offer
an update when that version changes. Never put a `version` into a plugin's own
`plugin.json`: Claude Code would let it silently override the marketplace
version.
