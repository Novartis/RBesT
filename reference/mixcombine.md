# Combine Mixture Distributions

Combining mixture distributions of the same class to form a new mixture.

## Usage

``` r
mixcombine(..., weight, rescale = TRUE)
```

## Arguments

- ...:

  arbitrary number of mixtures with same distributional class. Each
  component with values of mixture weight and model parameters.

- weight:

  relative weight for each component in new mixture distribution. The
  vector must be of the same length as input mixtures components. The
  default value gives equal weight to each component.

- rescale:

  boolean value indicates if the weights are rescaled to sum to 1.

## Value

A R-object with the new mixture distribution.

## Details

Combines mixtures of the same class of random variable to form a new
mixture distribution.

## See also

[`robustify()`](https://opensource.nibr.com/RBesT/reference/robustify.md)

Other mixdist:
[`mix`](https://opensource.nibr.com/RBesT/reference/mix.md),
[`mixbeta()`](https://opensource.nibr.com/RBesT/reference/mixbeta.md),
[`mixgamma()`](https://opensource.nibr.com/RBesT/reference/mixgamma.md),
[`mixjson`](https://opensource.nibr.com/RBesT/reference/mixjson.md),
[`mixmvnorm()`](https://opensource.nibr.com/RBesT/reference/mixmvnorm.md),
[`mixnorm()`](https://opensource.nibr.com/RBesT/reference/mixnorm.md),
[`mixplot`](https://opensource.nibr.com/RBesT/reference/mixplot.md)

## Examples

``` r
# beta with two informative components
bm <- mixbeta(inf = c(0.5, 10, 100), inf2 = c(0.5, 30, 80))

# robustified with mixcombine, i.e. a 10% uninformative part added
unif <- mixbeta(rob = c(1, 1, 1))
mixcombine(bm, unif, weight = c(9, 1))
#> Univariate beta mixture
#> Mixture Components:
#>   inf    inf2   rob   
#> w   0.45   0.45   0.10
#> a  10.00  30.00   1.00
#> b 100.00  80.00   1.00
```
