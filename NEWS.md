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
