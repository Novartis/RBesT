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
#> 1  0.207  0.207 0.0416 0.123 0.207 0.291
#> 2  0.203  0.203 0.0384 0.129 0.203 0.280
#> 3  0.221  0.218 0.0349 0.157 0.218 0.297
#> 4  0.242  0.238 0.0360 0.180 0.238 0.322
#> 5  0.200  0.200 0.0287 0.144 0.200 0.257
#> 6  0.186  0.184 0.0387 0.114 0.184 0.263
#> 7  0.226  0.222 0.0406 0.157 0.222 0.316
#> 8  0.175  0.175 0.0385 0.100 0.175 0.252
#> 9  0.228  0.223 0.0518 0.135 0.223 0.351
#> 10 0.182  0.182 0.0349 0.116 0.182 0.250
#> 11 0.235  0.234 0.0277 0.184 0.234 0.292

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
#> 1  0.218  0.213 0.0604 0.1130 0.213 0.357
#> 2  0.217  0.213 0.0582 0.1130 0.213 0.347
#> 3  0.217  0.212 0.0610 0.1130 0.212 0.359
#> 4  0.219  0.215 0.0610 0.1130 0.215 0.357
#> 5  0.219  0.214 0.0596 0.1160 0.214 0.364
#> 6  0.198  0.193 0.0603 0.0933 0.193 0.343
#> 7  0.198  0.193 0.0597 0.0972 0.193 0.337
#> 8  0.198  0.193 0.0625 0.0933 0.193 0.339
#> 9  0.219  0.214 0.0652 0.1070 0.214 0.369
#> 10 0.217  0.212 0.0643 0.1040 0.212 0.369
#> 11 0.217  0.211 0.0653 0.1070 0.211 0.381


# a summary function returns the results as matrix
summary(pred_cov)
#>         mean    median         sd      q2.5       q50     q97.5
#> 1  0.2065430 0.2068254 0.04164911 0.1231897 0.2068254 0.2906691
#> 2  0.2032181 0.2032065 0.03844876 0.1288859 0.2032065 0.2798751
#> 3  0.2206000 0.2182938 0.03486121 0.1567846 0.2182938 0.2969690
#> 4  0.2418891 0.2384320 0.03597107 0.1804723 0.2384320 0.3219211
#> 5  0.1997706 0.2000837 0.02866866 0.1439139 0.2000837 0.2565522
#> 6  0.1855770 0.1844776 0.03869977 0.1141460 0.1844776 0.2634910
#> 7  0.2261279 0.2220795 0.04056395 0.1568187 0.2220795 0.3161546
#> 8  0.1751656 0.1752809 0.03848595 0.1000778 0.1752809 0.2519715
#> 9  0.2276877 0.2229155 0.05180018 0.1349640 0.2229155 0.3506505
#> 10 0.1822711 0.1820796 0.03487839 0.1158106 0.1820796 0.2501153
#> 11 0.2354901 0.2342259 0.02767388 0.1842712 0.2342259 0.2922140

# obtain a prediction for new data with specific covariates
pred_new <- predict(map, data.frame(country = "CH", study = 12))
pred_new
#> Meta-Analytic-Predictive Prior Predictions
#> Scale: response 
#> 
#> Summary:
#>    mean median     sd  q2.5   q50 q97.5
#> 1 0.218  0.213 0.0611 0.113 0.213 0.366
## Recover user set sampling defaults
options(.user_mc_options)
```
