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
.user_mc_options <- options()

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

# posterior predictive distribution for each input data item (shrinkage estimates)
pred_cov <- predict(map)
pred_cov
#> Meta-Analytic-Predictive Prior Predictions
#> Scale: response 
#> 
#> Summary:
#>     mean median     sd  q2.5   q50 q97.5
#> 1  0.206  0.206 0.0417 0.121 0.206 0.293
#> 2  0.202  0.203 0.0393 0.123 0.203 0.280
#> 3  0.222  0.220 0.0350 0.159 0.220 0.297
#> 4  0.242  0.239 0.0354 0.180 0.239 0.319
#> 5  0.200  0.200 0.0283 0.144 0.200 0.256
#> 6  0.185  0.185 0.0390 0.110 0.185 0.264
#> 7  0.226  0.222 0.0405 0.157 0.222 0.314
#> 8  0.175  0.175 0.0391 0.101 0.175 0.251
#> 9  0.227  0.222 0.0494 0.138 0.222 0.339
#> 10 0.183  0.183 0.0341 0.117 0.183 0.248
#> 11 0.236  0.235 0.0271 0.185 0.235 0.294

# extract sample as matrix
samp <- as.matrix(pred_cov)

# predictive distribution for each input data item (if the input studies were new ones)
pred_cov_pred <- predict(map, trans_cov)
pred_cov_pred
#> Meta-Analytic-Predictive Prior Predictions
#> Scale: response 
#> 
#> Summary:
#>     mean median     sd   q2.5   q50 q97.5
#> 1  0.218  0.213 0.0600 0.1170 0.213 0.364
#> 2  0.217  0.213 0.0583 0.1100 0.213 0.349
#> 3  0.217  0.212 0.0596 0.1110 0.212 0.357
#> 4  0.220  0.213 0.0610 0.1170 0.213 0.364
#> 5  0.219  0.215 0.0599 0.1110 0.215 0.363
#> 6  0.201  0.195 0.0642 0.0919 0.195 0.349
#> 7  0.201  0.194 0.0628 0.0958 0.194 0.351
#> 8  0.200  0.195 0.0621 0.0922 0.195 0.341
#> 9  0.220  0.213 0.0671 0.1060 0.213 0.387
#> 10 0.218  0.213 0.0642 0.1050 0.213 0.369
#> 11 0.218  0.213 0.0629 0.1060 0.213 0.364


# a summary function returns the results as matrix
summary(pred_cov)
#>         mean    median         sd      q2.5       q50     q97.5
#> 1  0.2062488 0.2062227 0.04174574 0.1214406 0.2062227 0.2930378
#> 2  0.2020430 0.2031260 0.03926258 0.1228433 0.2031260 0.2802455
#> 3  0.2215977 0.2195194 0.03497005 0.1588652 0.2195194 0.2972463
#> 4  0.2423708 0.2392820 0.03543783 0.1801468 0.2392820 0.3189910
#> 5  0.1995879 0.1999883 0.02828960 0.1440262 0.1999883 0.2558488
#> 6  0.1853129 0.1848818 0.03901433 0.1101961 0.1848818 0.2642492
#> 7  0.2259685 0.2224310 0.04054221 0.1570393 0.2224310 0.3135865
#> 8  0.1748746 0.1747496 0.03906731 0.1011184 0.1747496 0.2511657
#> 9  0.2269397 0.2224007 0.04939518 0.1382642 0.2224007 0.3388369
#> 10 0.1828821 0.1825474 0.03411280 0.1167423 0.1825474 0.2478503
#> 11 0.2357563 0.2346353 0.02714248 0.1850796 0.2346353 0.2937389

# obtain a prediction for new data with specific covariates
pred_new <- predict(map, data.frame(country = "CH", study = 12))
pred_new
#> Meta-Analytic-Predictive Prior Predictions
#> Scale: response 
#> 
#> Summary:
#>    mean median     sd  q2.5   q50 q97.5
#> 1 0.218  0.213 0.0598 0.112 0.213 0.359
## Recover user set sampling defaults
options(.user_mc_options)
```
