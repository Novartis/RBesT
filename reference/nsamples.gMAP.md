# Return the number of posterior samples

Return the number of posterior samples

## Usage

``` r
# S3 method for class 'gMAP'
nsamples(object, ...)
```

## Arguments

- object:

  fitted model object

- ...:

  not used in this function

## Examples

``` r
.user_mc_options <- options()


set.seed(34563)
map_AS <- gMAP(cbind(r, n - r) ~ 1 | study,
  family = binomial,
  data = AS,
  tau.dist = "HalfNormal", tau.prior = 1,
  beta.prior = 2
)
#> Assuming default prior location   for beta: 0

nsamples(map_AS)
#> [1] 4000

## Recover user set sampling defaults
options(.user_mc_options)
```
