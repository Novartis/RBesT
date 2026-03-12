# Predictions from gMAP analyses

Produces a sample of the predictive distribution.

## Usage

``` r
# S3 method for class 'gMAP'
predict(
  object,
  newdata,
  type = c("response", "link"),
  probs = c(0.025, 0.5, 0.975),
  na.action = na.pass,
  thin,
  ...
)

# S3 method for class 'gMAPpred'
print(x, digits = 3, ...)

# S3 method for class 'gMAPpred'
summary(object, ...)

# S3 method for class 'gMAPpred'
as.matrix(x, ...)
```

## Arguments

- newdata:

  data.frame which must contain the same columns as input into the gMAP
  analysis. If left out, then a posterior prediction for the fitted data
  entries from the gMAP object is performed (shrinkage estimates).

- type:

  sets reported scale (`response` (default) or `link`).

- probs:

  defines quantiles to be reported.

- na.action:

  how to handle missings.

- thin:

  thinning applied is derived from the `gMAP` object.

- ...:

  ignored.

- x, object:

  gMAP analysis object for which predictions are performed

- digits:

  number of displayed significant digits.

## Details

Predictions are made using the \\\tau\\ prediction stratum of the gMAP
object. For details on the syntax, please refer to
[`predict.glm()`](https://rdrr.io/r/stats/predict.glm.html) and the
example below.

## See also

[`gMAP()`](https://opensource.nibr.com/RBesT/reference/gMAP.md),
[`predict.glm()`](https://rdrr.io/r/stats/predict.glm.html)

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 20x more warmup & iter in practice
.user_mc_options <- options(RBesT.MC.warmup=50, RBesT.MC.iter=100,
                            RBesT.MC.chains=2, RBesT.MC.thin=1)

# create a fake data set with a covariate
trans_cov <- transform(
  transplant,
  country = cut(1:11, c(0, 5, 8, Inf), c("CH", "US", "DE"))
)
set.seed(34246)
map <- gMAP(
  cbind(r, n - r) ~ 1 + country | study,
  data = trans_cov,
  tau.dist = "HalfNormal",
  tau.prior = 1,
  # Note on priors: we make the overall intercept weakly-informative
  # and the regression coefficients must have tighter sd as these are
  # deviations in the default contrast parametrization
  beta.prior = rbind(c(0, 2), c(0, 1), c(0, 1)),
  family = binomial,
  ## ensure fast example runtime
  thin = 1,
  chains = 1
)
#> Warning: The largest R-hat is 1.2, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Maximal Rhat > 1.1. Consider increasing RBesT.MC.warmup MCMC parameter.
#> Final MCMC sample equivalent to less than 1000 independent draws.
#> Please consider increasing the MCMC simulation size.

# posterior predictive distribution for each input data item (shrinkage estimates)
pred_cov <- predict(map)
pred_cov
#> Meta-Analytic-Predictive Prior Predictions
#> Scale: response 
#> 
#> Summary:
#>     mean     sd  2.5%   50% 97.5%
#> 1  0.217 0.0429 0.166 0.206 0.308
#> 2  0.195 0.0490 0.134 0.186 0.271
#> 3  0.222 0.0419 0.150 0.220 0.298
#> 4  0.254 0.0430 0.190 0.252 0.334
#> 5  0.203 0.0257 0.157 0.201 0.249
#> 6  0.174 0.0383 0.110 0.178 0.243
#> 7  0.224 0.0439 0.140 0.225 0.306
#> 8  0.175 0.0394 0.101 0.170 0.245
#> 9  0.230 0.0505 0.155 0.229 0.349
#> 10 0.178 0.0350 0.118 0.181 0.235
#> 11 0.237 0.0319 0.167 0.239 0.288

# extract sample as matrix
samp <- as.matrix(pred_cov)

# predictive distribution for each input data item (if the input studies were new ones)
pred_cov_pred <- predict(map, trans_cov)
pred_cov_pred
#> Meta-Analytic-Predictive Prior Predictions
#> Scale: response 
#> 
#> Summary:
#>     mean     sd   2.5%   50% 97.5%
#> 1  0.222 0.0618 0.1250 0.223 0.335
#> 2  0.229 0.0567 0.1460 0.229 0.351
#> 3  0.223 0.0635 0.1260 0.214 0.350
#> 4  0.215 0.0643 0.1180 0.203 0.359
#> 5  0.220 0.0671 0.0907 0.219 0.373
#> 6  0.191 0.0700 0.1100 0.182 0.302
#> 7  0.199 0.0783 0.1160 0.180 0.402
#> 8  0.188 0.0508 0.1130 0.184 0.325
#> 9  0.234 0.0878 0.1200 0.209 0.472
#> 10 0.230 0.0707 0.1010 0.223 0.351
#> 11 0.220 0.0602 0.1260 0.225 0.329


# a summary function returns the results as matrix
summary(pred_cov)
#>         mean         sd      2.5%       50%     97.5%
#> 1  0.2170247 0.04292502 0.1659672 0.2059520 0.3077934
#> 2  0.1948558 0.04896236 0.1340393 0.1860567 0.2713109
#> 3  0.2219566 0.04191868 0.1497590 0.2197390 0.2977980
#> 4  0.2535315 0.04295406 0.1898418 0.2518005 0.3337547
#> 5  0.2033785 0.02570938 0.1568488 0.2011043 0.2485448
#> 6  0.1739990 0.03833674 0.1097742 0.1781922 0.2429666
#> 7  0.2239574 0.04388595 0.1395464 0.2245863 0.3058693
#> 8  0.1746536 0.03936547 0.1007810 0.1701036 0.2447119
#> 9  0.2304102 0.05048495 0.1551377 0.2288809 0.3486094
#> 10 0.1778153 0.03502150 0.1180302 0.1812999 0.2346684
#> 11 0.2368112 0.03194805 0.1666561 0.2390441 0.2876846

# obtain a prediction for new data with specific covariates
pred_new <- predict(map, data.frame(country = "CH", study = 12))
pred_new
#> Meta-Analytic-Predictive Prior Predictions
#> Scale: response 
#> 
#> Summary:
#>    mean     sd  2.5%   50% 97.5%
#> 1 0.227 0.0656 0.111 0.218 0.348
## Recover user set sampling defaults
options(.user_mc_options)
```
