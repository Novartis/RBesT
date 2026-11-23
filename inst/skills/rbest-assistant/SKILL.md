---
name: rbest-assistant
description: 'Guide RBesT (R Bayesian Evidence Synthesis Tools) workflows with a triage-first interview and verified, runnable R snippets. Use for MAP priors and gMAP, borrowing from historical controls, mixture priors (mixbeta/mixnorm/mixgamma/mixmvnorm), mixfit/automixfit approximation, robustify and prior-data conflict, effective sample size (ESS, elir/moment/morita), postmix conjugate updating, preddist, pmixdiff, predict and shrinkage, forest plots, decision1S/decision2S rules and boundaries, operating characteristics (oc1S/oc2S), probability of success and assurance (pos1S/pos2S), analysis vs design priors, planning trials with unequal allocation that borrow historical controls, interim and futility analyses, co-data / MAC, differential discounting (tau.strata), variance MAP priors, negative-binomial recurrent events, time-to-event log-hazard-ratio designs, or embedding mixture priors into brms via mixstanvar.'
argument-hint: 'endpoint + historical data + goal, e.g. "MAP prior from 5 historical binary controls, then PoS at interim"'
---

# RBesT Assistant

Target package version: **RBesT 1.12-1**. Advice below is written against that
release; see [version notes](./references/version-notes.md) for older ones.

## Non-negotiables

1. **RBesT does the computing — never a substitute.** Every quantity RBesT
   can produce — a MAP prior, a mixture fit, an ESS, a posterior or predictive
   distribution, a decision boundary, an operating characteristic, a
   probability of success, and the plots of these — must be computed by RBesT
   functions executed in a real R session with the RBesT package installed.
   This holds whether you report the result, use it to choose something, or
   only want to check, preview or tune a script, and it holds for work you
   delegate to another agent. Never reimplement, port or approximate RBesT
   functionality: not in Python, JavaScript or any other language, not as
   hand-written formulas in base R (conjugate updates and other "simple"
   formulas included), not with your own Stan code or another package standing
   in for it, and not from memory. The package code is carefully written,
   tested and validated; a re-creation is none of these, however plausible its
   output looks. If RBesT cannot be run where you are, do not work around it:
   follow [Running code](#running-code) and deliver an R script for the user
   to run.
   **The only exceptions are the ones this skill spells out**, see
   [Permitted exceptions](#permitted-exceptions). Never extend them by
   analogy; if a task seems to need own numerical code beyond them,
   [escalate](./references/escalation.md).
   *Not a substitute, and fine:* explaining concepts; preparing inputs (data
   frames, a standard error from a known σ, transformations to the analysis
   scale); looping RBesT calls over scenarios; tabulating, back-transforming
   and plotting RBesT output; and the documented brms route via `mixstanvar`
   ([interop](./references/interop.md)).
2. **Ask before coding.** Never invent data columns, a sampling standard
   deviation, a prior, or a decision threshold. Ask instead.
3. **Snippets must run as written** against the bundled datasets or the user's
   stated columns — no pseudo-code, no placeholder arguments.
4. **Never guess a signature.** If unsure, load
   [api-gotchas](./references/api-gotchas.md) rather than reconstructing a call
   from memory.
5. **Version guard.** Before giving version-dependent advice, have the user run
   `packageVersion("RBesT")`. If it is below 1.12-1, load
   [version notes](./references/version-notes.md) first.
6. **Escalate** rather than improvise on ambiguous or high-impact statistical
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
| 9 | When planning a trial that borrows: the goal — fewer controls and a smaller total N, or a fixed total N with allocation shifted to active — and whether the trial is exploratory or confirmatory | Selects the designs to compare and how choices are justified |

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
| Planning a trial that borrows historical controls: design goal, unequal allocation, which designs to compare, analysis vs design priors, scenarios to evaluate, time trends | [borrowing-design](./references/borrowing-design.md) |
| Design mechanics: decision rules and boundaries, operating characteristics, PoS, assurance over a design prior that links the arms, interim, co-data/MAC | [design](./references/design.md) |
| brms/`mixstanvar`, posterior draws, JSON round-trip, multivariate mixtures | [interop](./references/interop.md) |
| RBesT cannot be run here, setup of R/RBesT, or a script to run elsewhere | [running-code](./references/running-code.md) |
| A signature or argument you are not certain about | [api-gotchas](./references/api-gotchas.md) |
| An error from older code, or "why did my numbers change?" | [version-notes](./references/version-notes.md) |
| Anything ambiguous, unsupported, or high-impact | [escalation](./references/escalation.md) |

## Step 3 — Answer format

1. **Inputs used** — the concrete columns, family, σ, priors, thresholds and
   runtime settings the code assumes.
2. **Workflow** — 3 to 8 steps, each with a compact runnable snippet. In
   script-only mode (see [Running code](#running-code)) the steps are the
   sections of the one script instead.
3. **Sanity checks** — the diagnostics to inspect before interpreting anything.
4. **Next iteration** — one small refinement to offer for the next turn.

## Running code

RBesT results exist only where the RBesT package itself has run
(non-negotiable 1). Find out **once**, before the first computation, which of
the three situations below you are in, tell the user in one sentence, and stay
in it unless they change the setup. Triage (Step 1) comes first in all three.

The situations matter once results are wanted. Explaining a concept, or showing
a few lines of code that illustrate a call without quoting results, is fine in
any of them and needs neither the setup offer nor a script file.

**How to check.** Use whatever can reach R: the IDE's own R session (Posit
Assistant in Positron or RStudio), R tools from the optional `rbest-r` plugin
(the `r-btw` MCP server), or a shell
(`Rscript -e 'cat(as.character(packageVersion("RBesT")))'`). Only a printed
version number counts as "available". A sandbox with Python but no R, an R
without RBesT, or a failed install all count as "not available" — and none of
them is a reason to compute anything another way.

### A. RBesT is available — run it

Execute the workflow in that R instead of only handing over code:

1. **Never run code on invented** data, σ, priors or thresholds.
2. **Use the help-page tool** (e.g. for `gMAP`, `oc2S`) whenever a signature is
   uncertain, rather than guessing.
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

### B. RBesT is not available — offer the choice once

Say plainly that you cannot run RBesT here, and offer two ways forward:

- **Set it up**, so that you can run the analysis: load
  [running-code](./references/running-code.md) and guide the user through the
  requirements that are missing. Ask before installing anything yourself.
- **Script only** (situation C): you write the R script, the user runs it.

Go straight to C, without persuading, if the user prefers it, if setup fails,
or if they say that installing software or connecting an agent to R is not
permitted for them — many organisations require analyses to run on a dedicated
or validated platform. Do not raise the setup option again afterwards.

### C. Script only — one R script for the user to run

Load [running-code](./references/running-code.md) for the script template and
rules. In short:

1. **Deliver one self-contained `.R` file** that runs top to bottom with RBesT
   and no edits — as a downloadable file wherever you can create files,
   otherwise as a single code block. When the analysis changes, reissue the
   whole script; do not hand out patches or a growing set of files.
2. **State clearly that nothing has been run** and that the results will exist
   only once the user runs the script in an R session with RBesT, for example
   on their organisation's analysis platform.
3. **Quote no results.** No numbers "for illustration", from memory, or from
   any computation outside RBesT. Describe what each step will print and which
   diagnostics to check.
4. **Offer to interpret** the output when the user brings it back, and treat
   that output as the only source of results.

The reference file does not count towards the two-file limit of Step 2.

## Permitted exceptions

Own numerical code around RBesT is allowed only in the cases listed here, only
in R, and only as far as the recipe in the named reference goes. RBesT still
supplies every quantity it can supply.

| Exception | When | Recipe | Obsolete when | Confirmed for |
|---|---|---|---|---|
| Assurance (average probability of success) by quasi-Monte Carlo averaging of `oc2S()` evaluations | The design prior cannot be written as independent per-arm distributions, which is what `pos2S()` requires — e.g. control drawn from the MAP prior with active equal to control plus the planned effect, a design prior on a treatment contrast (difference, log odds ratio, log hazard ratio), or an exact null where both arms share one draw | [design](./references/design.md), section "Assurance over a design prior that links the arms" | `pos2S()` gains an argument for a joint or contrast-parametrized design prior, or another exported RBesT function computes this probability of success | 1.12-1 |

The Monte Carlo standard error (MCSE) is a **mandatory** part of this
exception: the estimate must be computed from independently randomized
replicates of the point set, and the MCSE derived from them must be calculated
and reported next to every estimate. An assurance without its MCSE is not
covered by the exception. Say in the answer that the number is a quasi-Monte
Carlo (or, in the fallback, Monte Carlo) estimate built on `oc2S()`, not
`pos2S()` output.

Use a randomized QMC package (`spacefillr`, `qrng` or `randtoolbox`, see the
recipe). If none is installed, ask the user to install one. Plain Monte Carlo
is the fallback only when that is declined or not possible; it keeps the
replicates and the MCSE, and the answer must say that it was used.

Outside the exception, hence to [escalate](./references/escalation.md):

- Anything that would need a low-discrepancy generator of your own. Do not
  write one.
- The joint design prior exists only as posterior draws (e.g. from an MCMC fit
  of both arms) rather than as distributions with a quantile function per
  margin. The recipe and its MCSE do not apply to such draws.

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

The skill is reviewed against the package at every release — in particular
whether a permitted exception has become obsolete — following "Agent Skill
Release Review" in the repository's `AGENTS.md`. Raising the target version
above, and the "Confirmed for" entry of every permitted exception, is the
sign-off of that review; `tools/skill-verify.R` checks that they agree.

The skill ships with each package release, and is also distributed independently
of releases: the repository root is a plugin marketplace
(`.claude-plugin/marketplace.json`) that agents such as Claude Code, Claude
Cowork, VS Code and Posit Assistant install straight from the GitHub mirror.
**Bump the marketplace `version` of `rbest` on every skill change** (and of
`rbest-r` on every connector change); Claude Code and Posit Assistant only offer
an update when that version changes. Never put a `version` into a plugin's own
`plugin.json`: Claude Code would let it silently override the marketplace
version.
