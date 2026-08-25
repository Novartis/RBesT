# RBesT 1.12-0 - unreleased

## Enhancements

* `gMAP` now samples the group random effects in a sum-to-zero
  parametrization whenever the model has an intercept, normal random
  effects and a single `tau` stratum. The data-free common shift shared
  by the intercept and the random effects is marginalized analytically
  using an orthonormal basis of the sum-to-zero subspace, which removes
  the ridge between `beta[1]` and `mean(eps)` that made the posterior
  hard to explore. The model is mathematically unchanged and all outputs
  keep their previous meaning; only the sampling geometry differs.
  Student-t random effects, several `tau` strata and intercept-free
  models continue to use the previous parametrization. Set
  `options(RBesT.MC.s2z = FALSE)` to return to the previous
  parametrization.
* Lower the default `adapt_delta` for `gMAP` from `0.99` to `0.95`
  (`options(RBesT.MC.control)`). The sum-to-zero parametrization no
  longer needs the very conservative target acceptance rate that the old
  geometry required, and the lower target draws the same effective
  sample size from substantially fewer gradient evaluations. Should
  divergent transitions still occur, the reported warning now recommends
  raising `adapt_delta` to `0.99`. Setting
  `options(RBesT.MC.s2z = FALSE)` also restores the previous
  `adapt_delta` default of `0.99`, since the legacy geometry needs the
  more conservative target acceptance rate.

# RBesT 1.11-0 - August 3rd, 2026

## Enhancements

* Decision boundary estimation is now substantially faster. The fixed
  uniform grid sweep in `decision2S_boundary` / `decision1S_boundary`
  (and hence `oc1S`, `oc2S`, `pos1S`, `pos2S`) is replaced by an adaptive
  scheme that only solves where the boundary bends (normal) or fills
  provably constant runs of the monotone integer boundary without further
  root searches (binomial, Poisson). Boundary construction is roughly
  `10x` faster for normal endpoints and up to orders of magnitude faster
  for near-flat discrete boundaries (e.g. large-`n` Poisson); normal
  `oc2S`/`pos2S` run about `2-3x` faster. Boundary values are unchanged
  (identical for discrete, within `~1e-5` for normal). Set
  `options(RBesT.decision2S_boundary = "grid")` to restore the legacy
  sweep.
* Make the normal `decision2S_boundary` interpolant consistent with the
  adaptive point placement. Curved and mixture boundaries are now
  interpolated with a shape-preserving monotone Hermite spline
  (`splinefun(method = "monoH.FC")`) instead of a global cubic, which
  matches the adaptive tracer's curvature-graded points and avoids
  overshoot. For the provably-linear case (single-component priors with
  constant `sigma`, including multi-criteria and `n2 = 0` decisions) the
  boundary is now built as an exact straight line from two solves rather
  than a grid, making these boundaries exact (previously accurate to
  `~1e-6`). Operating characteristics / probability of success are
  unchanged to within tolerance.
* Lower the default decision boundary discretization density `Ngrid`
  from `10` to `5` in `oc2S`, `pos2S`, and `decision2S_boundary` for the
  normal case. `Ngrid` now sets the finest resolution of the adaptive
  scheme and the grid density of the `"grid"` fallback; for the grid path
  a coarser grid roughly halves the boundary cost. Either way the
  resulting operating characteristics / probability of success change by
  at most about `1e-5` across a range of mixture and exponential-family
  settings. Pass `Ngrid = 10` to restore the previous behaviour.
* Consolidate all package references into a single BibTeX file at
  `inst/REFERENCES.bib`. Documentation now cites it via the `Rdpack`
  `\insertRef` macro (roxygen `@references`), and the vignettes and
  articles cite the same file through the R Markdown `bibliography:`
  field, removing duplicated per-file reference lists.

## Bugfixes

* Fix the simulation based calibration (SBC) code under `inst/sbc`,
  which was shipped broken in 1.10-0: the `gMAP` posterior storage
  refactor removed the raw `rstan` `stanfit` `fit` slot, causing
  `run_sbc_case` to fail with `no applicable method for '@' applied to
  an object of class "NULL"`. The SBC helper `fit_rbest` now reads
  sampler diagnostics and posterior draws from the stored `posterior`
  draws objects instead of the removed `stanfit`.

# RBesT 1.10-0 - June 30th, 2026

## Enhancements

* Add a `family` argument to `oc1S`, `oc2S`, `pos1S`, `pos2S`,
  `decision1S_boundary`, and `decision2S_boundary` for the normal
  case. Under the normal approximation this allows operating
  characteristics, probability of success, and decision boundaries to
  be evaluated for arbitrary exponential family endpoints.
* Add a new article (vignette) demonstrating a futility interim
  analysis for a negative binomial endpoint, including grid-based
  operating characteristics for two-stage designs, MAP-prior borrowing
  on control at the interim, and conditional vs. predictive power.
* Use Gaussian quadrature (GQ) for numerical integration against
  mixture densities (Gauss-Hermite for `normMix`, Gauss-Jacobi for
  `betaMix`, Gauss-Laguerre for `gammaMix`), replacing adaptive
  Gauss-Kronrod integration on identity-link paths in `oc2S`, `pos2S`,
  `dmixdiff`, `pmixdiff`, and `ess` for deterministic and robust
  results. The node count is refined automatically (starting from
  `RBesT.GQ_nodes`) until successive estimates agree within
  `max(RBesT.GQ_abs_tol, RBesT.GQ_rel_tol * |I|)`, falling through to
  adaptive integration if the node cap is reached. Non-identity link
  cases and betaMix ESS continue to use adaptive integration, which can
  also be restored globally with
  `options(RBesT.integrate_method = "adaptive")`.
* Use modern posterior functions to produce posterior summaries. Note
  that this can change the numeric values wrt to previous versions
  (e.g. a modernized Rhat definition is used).
* Refactor `gMAP` posterior sample storage to use `posterior` draws
  objects plus aligned sampler diagnostics instead of relying on a raw
  `rstan` `stanfit`; update summaries, predictions, plotting, and
  posterior accessors accordingly. The previously stored `fit` slot
  has been removed.
* Enable default test coverage for `gMAP` related functions using
  compact approximate posterior draw fixtures, so these tests can run
  routinely without full MCMC fixture generation.

## Bugfixes

* Make `read_mix_json` robust against `jsonlite` 2.0.0 simplification
  behavior by using deterministic JSON parsing, fixing reads of normal
  mixtures with three or more components.

# RBesT 1.9-0 - March 13th, 2026

## Enhancements

* Reformat R sources using `Air`.
* Extend normal outcome functions to two-sided decisions: `decision1S`
  and `decision2S` now allow `lower.tail` to have as many elements as
  `pc` to allow the specification of two-sided decision boundaries
  (mixed `lower.tail` elements, i.e. some specifying "lower" and some
  "upper") to capture intermediate result scenarios.
* Introduce global option `RBesT.MC.save_warmup`, which is set to
  `FALSE` causing to drop warmup samples from the `gMAP` fit
  object. Note that the pre-1.9.0 behavior was to store the warmup
  samples (but these were not used). To restore the original behavior
  set `options(RBesT.MC.save_warmup=TRUE)` once globally.
* Add experimental `as_draws*` functions allowing to extract samples
  in the posterior draws format from `gMAP` objects.
  
## Bugfixes

* Avoid use of `size` argument to `ggplot2` routines, which have been
  deprecated.

# RBesT 1.8-2 - April 25th, 2025

## Enhancements

* Add experimental `write_mix_json` and `read_mix_json` functions
  which write and read mixture objects as JSON to/from files.
* Add `"msr"` parametrization to `mixmvnorm` which allows to consturct
  multi-variate normal mixtures in a parametrization using the mean
  vector, standard deviations and the correlations. This functionality
  is based on the new `msr2mvnorm` utility function.
* Allow standard deviations of 0 for components of a multi-variate
  normal mixture.

# RBesT 1.8-1 - January 20th, 2025

## Bugfixes

* Fixed an issue with the `ess` function for beta and gamma mixtures
  when used inside of `lapply` or `sapply`.

# RBesT 1.8-0 - January 8th, 2025

## Enhancements

* Enable ESS calculation for normal mixture densities when used in the
  context of a standard one-parameter exponential family through the
  new `family` argument. For example, this can be used to calculate
  the ESS of a normal mixture density representing a logit transformed
  response scale.
* Reformat R sources using `styler`.

## Bugfixes

* Correct boundary behavior of `BinaryExactCI` function whenever no
  responses or no non-responses are observed. Fixes issue #21.
* Stabilize internal beta mixture information function, which corrects
  unstable ESS ELIR computations. Addresses issue #22.

# RBesT 1.7-4 - November 21st, 2024

## Enhancements

* Added for `mixstanvar` automatic generation of distribution
  functions allowing truncated priors in `brms` with mixtures
* Slight speed increase for Stan model by more efficient construction
  of likelihood. Also reducing object size of `gMAP` objects by
  avoiding to store redundant variables.
* Upgraded `testthat` edition to version 3.

## Bugfixes

* Fix issue #18 of ESS ELIR aborting whenever one mixture component
  has zero weight.
* Fixed a rare issue when estimating ESS ELIR for beta mixtures. The
  calculation aborted due to boundary issues which are now avoided.
* Avoid over and underflow of beta-binomial distribution.
* Replace calls to deprecated `ggplot2::qplot` with respective calls
  to `ggplot2::ggplot`.
* Use new `array` notation for `mixstanvar` generated Stan code to
  make generated Stan programs compatible with Stan >= 2.33

# RBesT 1.7-3 - January 2nd, 2024

## Enhancements

* Updated Stan model file syntax to use new array syntax as required
  by Stan >=2.33. This upgrades the minimal Stan version to 2.26.
* Moved most vignettes to be articles to decrease size of R
  packages. Articles are available on [the package
  homepage](https://opensource.nibr.com/RBesT/articles/).

## Bugfixes

* Added `pngquant` to system requirements of package as requested by
  CRAN.

# RBesT 1.7-2 - August 21st, 2023

## Enhancements

* Add unit tests for `mixstanvar` function and disabled integration
  tests on CI/CD systems as these require significant resources due to
  the need to compile Stan models. The unit tests do ensure that
  things are correct from the `RBesT` side of things.
* Use new scheme of `roxygen2` to document package documentation.

## Bugfixes

* Make EM for normal and multivariate normal mixtures deterministic in
  that the resulting density does not depend on the random seed any
  longer, which is the intended behavior. Note that this change may
  cause that estimated normal / MVN mixtures will be numerically
  different as compared to the previous version. The estimated density
  itself is (statistically) still the very same.

# RBesT 1.7-1 - August 8th, 2023

## Enhancements

* Allow use of mixture priors in `brms` models with the new adapter
  function `mixstanvar`.
* Allow use of named dimensions for multivariate normal mixtures.
* Add *experimental* diagnostic plots for mixture multivariate normal
  EM fits. These are subject to changes in future release and user
  feedback is very welcome.
* Compress png plots of vignettes saving ~40% in file size of the
  package sources.

## Bugfixes

* Fix issue when plotting EM diagnostic debugging plots for normal
  mixtures.
* Comply with newer ggplot2 standards to not use `aes_string`.

# RBesT 1.7-0 - July 19th, 2023

## Enhancements

* Implementation multivariate normal mixture support in a first
  version. This includes density evaluation, basic summary functions
  and multivariate normal EM fitting. Support is not yet as complete
  as for other densites, but will be expanded in upcoming releases.
* Change the default for the option `constrain_gt1` of the EM for beta
  mixtures to `TRUE`. This will be default constrain the fitted
  parameters of each beta mixture component to be greater than unity
  as required for finite ESS elir calculations.

# RBesT 1.6-7 - June 26th, 2023

## Bug fixes

* resolve compilation from source issues on some platforms triggered
  by changes in `rstantools`
* correct documentation on difference distribution and improve PoS
  with co-data vignette

# RBesT 1.6-6 - March 3rd, 2023

## Bug fixes

* ensure C++17 compatiblity per CRAN (triggers an issue with
  clang 16)
* fix links in README.md to now link to new public pkgdown web-site

# RBesT 1.6-5 - February 8th, 2023

## Enhancements

* upon package load `RBesT` will now report the date of the release
  and the respective git commit hash used to create the sources of the
  package.

## Bug fixes

* ensure that `predict` for new studies will sample the study specifc
  random effect per iteration only once. This is important for MAP
  priors with covariate effects (which are pooled over the studies).

# RBesT 1.6-4 - August 5th, 2022

## Enhancements

* use clustermq inplace of batchtools for SBC runs. Also use L'Ecuyer
  CMG random number gen during SBC runs.
* expand introductory vignette with plot of ESS vs robust weight
* start using `matrixStats` to speed up EM algorithms & OCs
* avoid warning whenever 2S normal design expanded their domain when
  called repeatedly
* add `mixdist` plot when plotting a mixture resulting from *mixfit call
* add warning message when printing a `gMAP` analysis for which the
  Stan sampler had issues due to divergent transitions or non-convergence
  as indicated by large Rhat

## Bug fixes

* mixture density evaluations (`dmix`) with a defined link function
  did not evaluate correctly, which was visible when plotting mixtures
  (with more than one component) with link functions

# RBesT 1.6-3 - November 23rd, 2021

* update references to JSS publication

# RBesT 1.6-2 - September 3rd, 2021

* link against RcppParallel to align with new Stan requirements
* address CRAN comments

# RBesT 1.6-1 - May 28th, 2020

* stabilize elir ESS integration by integrating per mixture component
* comply with forthcoming and stricter stanc3 Stan transpiler
* address some warnings from ggplot2 3.3.0

# RBesT 1.6-0 - March 27th, 2020

* fix CI system issues
* fix issues with normal decision2S_boundary when boundaries are grown
* add demo for 2S OC simulation code for time-to-event endpoint with
  constant hazard assumption
* drop tidyverse dependency
* expand SBC checks to include group specifc estimates
* stop setting the ggplot2 default theme when loading package. All
  plots now use the bayesplot theme which can be modified with
  bayesplot_theme_* functions. See ?bayesplot::bayesplot_theme_get.
* correct transformation issue in MAP for variances vignette - thanks
  to Ping Chen
* allow for constrained fitting of beta mixtures which have a & b
  parameters greater than 1. This is the new default behavior which
  the function will inform about. The informational message will be
  removed in a future release.
* introduced new mixecdf plot as diagnostic for EM fits

# RBesT 1.5-4 - October 22nd, 2019

* Now really fix n2=0 case for 2S design functions for indirect
  comparisons
* Update package structure to new rstantools 2.0 system

# RBesT 1.5-3 - August 28th, 2019

* Fix vignette MAP for variances (missing definition)

# RBesT 1.5-2 - August 28th, 2019

* Speedup example run time
* Avoid use of cat in functions and use message instead
* Replace dontrun in examples with donttest
* Require NOT_CRAN=true for tests to run

# RBesT 1.5-1 - August 28th, 2019

* Work around compiler warning with clang on fedora platform

# RBesT 1.5-0 - August 15th, 2019

* Fix indirect comparisons to work with normal/Poisson/binomial
  (inexact) to allow for n2=0 in oc2S calls.
* Make mixture quantile function more robust to work with very flat
  mixture priors.
* Align ESS Morita calculations with Neuenschwander B. et al.,
  _pre-print_ 2019; arXiv:1907.04185

# RBesT 1.4-0 - July 27th, 2019

* Introduce elir ESS method as new default for ESS
* Allow to sample prior predictive with gMAP (argument prior_PD)
* Switch internally to ab parametrized version of EM beta algorithm

# RBesT 1.3-8 - April 3nd, 2019

* Use Simulation Based Calibration for gMAP model qualification
* Improve covariate handling (naming of data items)
* Speedup Stan model by avoiding matrix-vector products with many zeros
* Fix index issue with differential discounting when used with covariates
* Make initialization of EM algorithms more robust
* Avoid special build hacks on MacOS

# RBesT 1.3-7 - November 16th, 2018

* Address issue for build process on MacOS.

# RBesT 1.3-6 - November 14th, 2018

* Re-create vignettes with proper MCMC sampling.
* Automate R package build process using CI/CD.

# RBesT 1.3-5 - November 13th, 2018

* Corrected 1.3-4 release notes to include MAP for variances vignette
* Make build process more robust (updated src/Makevars{.win})
* Added probability of success at interim vignette
* Added probability of success with co-data vignette

# RBesT 1.3-4 - October 16th, 2018

* Make package work with rstan 2.18.1.
* Revert BetaBinomial implementation back to R functions.
* Bugfix for decision1S_boundary for normal case for extreme parameter configurations (fixes pos1S & oc1S as well).
* Bugfix for mixcombine and plot with normal mixtures without a sigma being defined.
* Bugfix for repeated calls to decision2S_boundary for normal endpoint (fixes pos2S & oc2S as well).
* Avoid use of deprecated bayesplot function arguments whenever divergencies occured.
* (corrected) Added MAP for variances vignette

# RBesT 1.3-3 - February 2nd, 2018

* Change numerical equality testing to use expect_equal (which uses
  all.equal internally accounting for machine specifc tolerances) to
  pass tests for no long double case. Numerical tolerances are
  reverted back to 1.3-1 settings.

# RBesT 1.3-2 - January 25th, 2018

* Adjust numerical tolerances to pass tests for no long double case

# RBesT 1.3-1 - December 21st, 2017

* Add Trustees of Columbia copyright for respective files in
  DESCRIPTION

# RBesT 1.3-0 - December 21st, 2017

* Added probability of success calculation for 1+2 sample case.
* Added decision1+2S_boundary functions (and deprecated use of y2
  argument of oc functions)
* Added RBesT.integrate_args option for greater control over density
  integrations.
* Correct cumulative predictive of beta mixtures to return 0/1 for
  out-of-range values (instead of leaving those out).
* Deprecated functions oc1+2Sdecision which are replaced by
  decision1+2S.

# RBesT 1.2-3 - August 21st, 2017

* Fix plotting procedures to work with bayesplot 1.3.0

# RBesT 1.2-2 - July 12th, 2017

* Further speedup example runtimes.

# RBesT 1.2-1 - July 12th, 2017

* Compactify reference PDF manual.
* Introduce sampling arguments to gMAP.
* Shorten runtime of examples.

# RBesT 1.2-0 - July 3rd, 2017

* First CRAN release.
* Update of documentation.

# RBesT 1.1-0 - May 15th, 2017

* Redesign of reference scale handling for normal case.
* Enable standard error as sufficient statistic in \code{postmix} function.
* Introduced plotting options.
* Increased adapt_delta default and set stepsize+max_treedepth default.
* Added RBesT.MC.\{ncp, init, rescale\} option.
* Corrections for Poisson OC.
* pmixdiff function now integrates over full support.
* Added \code{link} argument to \code{oc2Sdecision} which enables designs based on log-odds decisions or relative risks.
* New graphical vignette and new forest plot function.
* Use \pkg{bayesplot} as standard plotting package.

# RBesT 1.0-0 - March 10th, 2016

* Stabilize integration in pmixdiff for beta mixtures by logit transform.
* Set default adapt_delta to 0.975.
* Made RBesT compatible with ggplot2 2.0.
* Allowed n2=0 in \code{oc2S} function.

# RBesT 0.9-2 - Oct 28th 2015

* Corrected Poisson stratified estimates.
* Added warning on divergent transitions.
* Added the crohn dataset.

# RBesT 0.9-1 - Sept 3rd 2015

* Minor typo fixes.

# RBesT 0.9-0 - Sept 1st 2015

* First release.
