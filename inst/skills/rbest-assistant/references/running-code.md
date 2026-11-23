# Running Code: Setup and Script-Only Delivery

Use when RBesT cannot be run where the agent is working, when the user wants
help getting R and RBesT in place, or when the deliverable is a script to be run
elsewhere.

The rule behind all of it: **results come from the RBesT package or not at
all.** If RBesT is not available, the choice is between making it available and
handing over a script. Computing the answer some other way is not on the list.

## What counts as a substitute

All of these are ruled out, including "just to give a first impression":

- porting RBesT functionality to Python, JavaScript, Julia or anything else the
  sandbox happens to offer;
- re-deriving the methods in base R — a hand-written EM for the mixture fit,
  hand-coded ESS, conjugate updates, decision boundaries or operating
  characteristics;
- writing your own Stan, JAGS or other model in place of `gMAP`;
- using another meta-analysis or Bayesian package as a stand-in;
- copying individual source files out of the RBesT repository and sourcing
  them instead of installing the package;
- quoting results from memory, from documentation examples presented as if they
  were the user's results, or from "typical" values.

The skill permits a small number of exceptions and names each one in
`SKILL.md` under "Permitted exceptions" — at present only the quasi-Monte Carlo
assurance built on `oc2S()` in [design](./design.md). They run in R on top of
RBesT, so they need RBesT just the same, and they are not a licence for
anything similar.

Installing the genuine package (from CRAN, or a development version from the
`Novartis/RBesT` GitHub repository when the user asks for one) is always fine —
that *is* RBesT.

## Finding out where you stand

Check with whatever reaches R, and trust only a printed version number:

```r
cat(as.character(packageVersion("RBesT")))
```

From a shell: `Rscript -e 'cat(as.character(packageVersion("RBesT")))'`.

Keep two places apart. An agent's sandbox (for example the virtual machine that
Claude Cowork works in) is not the user's computer: R found or missing there
says nothing about the user's machine, and the reverse. R tools that come from
an MCP server such as `r-btw` run in the user's R installation (a fresh R
process); the IDE's assistant (Posit Assistant) runs in the user's live R
session.

| Finding | Situation |
|---|---|
| A version number is printed | RBesT is available: run the analysis there |
| R answers, but RBesT is not installed | Offer setup (install the package) or script only |
| No R, or no tool that can run R at all | Offer setup (below) or script only |
| The user says running R through an agent is not permitted | Script only, no further discussion |

## Setup guidance

Give only the part the user is missing, and ask before installing anything
yourself. Stop and switch to script-only as soon as the user declines, a step
fails, or policy forbids it.

**R and RBesT on the user's computer**

Prefer a ready-built binary of RBesT over a build from source. A source build
compiles RBesT's Stan models together with large C++ dependencies: it needs a
C++ toolchain, a lot of memory and tens of minutes, and it is the usual point
of failure.

1. R from <https://cran.r-project.org> (RBesT needs R ≥ 3.5.0; use a current
   release, since binaries are built for current R versions only).
2. Install RBesT as a binary, preferably from Posit Public Package Manager
   (<https://packagemanager.posit.co>), which serves binaries for Windows,
   macOS and the common Linux distributions:

   - *Windows and macOS:*
     `install.packages("RBesT", repos = "https://packagemanager.posit.co/cran/latest")`.
     CRAN itself also provides Windows and macOS binaries, so a plain
     `install.packages("RBesT")` works as well.
   - *Linux:* plain CRAN installs from source. Binaries need the repository
     URL of the user's distribution, of the form
     `https://packagemanager.posit.co/cran/__linux__/<distribution>/latest`.
     Do not guess it: have the user take the URL, and the few lines of R
     configuration shown with it, from the setup page
     <https://packagemanager.posit.co/client/#/repos/cran/setup>.
   - *An organisation's own Posit Package Manager or internal repository*
     takes precedence over all of the above; use what the user's R is already
     configured for.

   Check that a binary was used: the install output must not show a
   compilation of RBesT. If it starts compiling, stop it and fix the
   repository setting rather than waiting.
3. Check: `library(RBesT); packageVersion("RBesT")`.

**Letting the agent use that R**

- *Positron or RStudio with Posit Assistant:* nothing more to do; it runs R in
  the live session.
- *Coding agents with a shell on the user's machine* (Claude Code, Copilot
  agent mode, Codex): `Rscript` must be on the `PATH`.
- *Claude Cowork and Claude Desktop:* their sandbox cannot reach the user's R.
  Install the optional `rbest-r` plugin from the same marketplace as this skill
  (`Novartis/RBesT`), plus `install.packages(c("btw", "ragg"))` in R, then
  restart the app. It starts btw's MCP server on macOS and Linux; for Windows
  and manual setups see <https://posit-dev.github.io/btw/reference/mcp.html>.
  The connector runs R code on the user's machine without a sandbox, and
  organisations may disable local MCP servers.

**Installing into the agent's own sandbox**

Only if R is already there and the user agrees: try one install of RBesT, from
a binary repository as described above wherever one fits the sandbox; say
beforehand that a source build takes tens of minutes and may run out of
memory. If there is no R, no network access to a package repository, no
compiler, or the install fails, stop there. Do not install R itself into a
sandbox: that is a larger intervention than this skill covers, and a missing R
means script-only or R on the user's side. Do not assemble a partial installation, and do not fall back on anything
listed under "What counts as a substitute".

**When setup is not possible or not allowed**

Accept it at once. Many organisations do not allow an agent to drive an R
session and require analyses to run on a dedicated or validated computing
platform. The script-only route is made for exactly this: the agent writes, the
user runs.

## Script-only delivery

The deliverable is **one R script** that the user downloads and runs in an R
session with RBesT.

**Rules**

1. One file, self-contained, runs top to bottom without edits — interactively
   or with `Rscript`. Name it after the task, e.g. `rbest_map_prior_binary.R`.
   Create it as a file the user can download whenever you can write files;
   otherwise give it as a single code block.
2. Triage first, as always. Every blocking input is settled before the script
   is written; nothing in it is a placeholder.
3. The data are part of the script: either written out as a `data.frame()` from
   the values the user gave, or read from the exact file path and columns they
   named. Bundled datasets (`AS`, `crohn`, …) only for illustration, and say so
   in the header.
4. Because you cannot test the script, be stricter than usual about signatures:
   take every call from this skill's reference files or
   [api-gotchas](./api-gotchas.md), not from memory. If an R without RBesT is
   at hand, you may check that the file parses (`parse(file = ...)`) and say
   that this is all that was checked. Do not "test" the script against a mock
   or a re-creation of RBesT functions.
5. The script stops with a clear message when RBesT is missing, sets
   `set.seed()`, prints every result under a label, prints the diagnostics that
   must be read before interpreting, and ends with `sessionInfo()`.
6. It only uses RBesT and the packages that come with R, unless the user asked
   for more (e.g. brms) or a recipe of this skill needs one (a QMC package
   for the assurance recipe — ask which one the user's platform has); check for such a package at
   the top, like for RBesT. It installs nothing and writes no files other than those announced in its
   header.
7. In the answer: say that nothing has been run, give no numbers, explain how
   to run the script and what to look at first, and offer to interpret the
   output if the user pastes or uploads it.
8. When the analysis changes, reissue the complete script under the same name.
   One current file, not a series of patches.

**Template.** Keep the frame; replace the numbered sections by the workflow
that triage arrived at.

```r
## ---------------------------------------------------------------------------
## RBesT analysis script: MAP prior for a binary endpoint (illustration)
##
## Written by an AI agent with the rbest-assistant skill. NOT EXECUTED by the
## agent: results exist only once you run this script yourself. Review code,
## inputs and results before use.
##
## Inputs and assumptions
##   data        : bundled dataset AS (8 studies; columns study, n, r)
##   model       : binomial, logit link, exchangeable studies
##   priors      : tau ~ HalfNormal(0, 1), intercept ~ Normal(0, 2^2)
##   robustify   : weight 0.2 on a vague component with mean 0.5
##   new data    : r = 10 responders out of n = 30
## Written for  : RBesT 1.12-1
## How to run   : source("rbest_map_prior_binary.R", echo = TRUE)
##                or: Rscript rbest_map_prior_binary.R > rbest_output.txt 2>&1
## Files written: none (under Rscript, plots go to Rplots.pdf)
## ---------------------------------------------------------------------------

if (!requireNamespace("RBesT", quietly = TRUE)) {
  stop("This script needs the RBesT package: install.packages(\"RBesT\")")
}
library(RBesT)
cat("RBesT version:", as.character(packageVersion("RBesT")), "\n")
if (packageVersion("RBesT") < "1.12.1") {
  warning("Script was written for RBesT 1.12-1; results may differ.")
}

set.seed(34563)

## 1. Data ------------------------------------------------------------------
hist_data <- AS
print(hist_data)

## 2. MAP prior (Stan MCMC; takes a moment) -----------------------------------
map_mcmc <- gMAP(cbind(r, n - r) ~ 1 | study,
  family = binomial, data = hist_data,
  tau.dist = "HalfNormal", tau.prior = 1,
  beta.prior = 2
)
cat("\n== MAP fit: check Rhat, divergences and tau before going on ==\n")
print(map_mcmc)
print(plot(map_mcmc)$forest_model)

## 3. Parametric mixture approximation ----------------------------------------
map_mix <- automixfit(map_mcmc)
cat("\n== Mixture approximation of the MAP prior ==\n")
print(map_mix)
print(plot(map_mix)$mix)

## 4. Robustification and effective sample size -------------------------------
map_robust <- robustify(map_mix, weight = 0.2, mean = 0.5)
cat("\n== Robust MAP prior ==\n")
print(map_robust)
cat("\n== ESS (elir): MAP, robust MAP ==\n")
print(c(map = ess(map_mix), robust = ess(map_robust)))

## 5. Posterior with the new data ---------------------------------------------
post <- postmix(map_robust, r = 10, n = 30)
cat("\n== Posterior ==\n")
print(summary(post))

## Session ------------------------------------------------------------------
cat("\n== Session ==\n")
print(sessionInfo())
```

**What to say alongside the script**

- "I could not run RBesT here, so this script has not been executed and I am
  not reporting any results."
- How to run it (the two lines from the header) and roughly how long the slow
  steps take.
- The diagnostics to read first, in the order the script prints them.
- "Paste or upload the output and I will go through it with you."
