# Beta Mixture Density

The Beta mixture density and auxilary functions.

## Usage

``` r
mixbeta(..., param = c("ab", "ms", "mn"))

ms2beta(m, s, drop = TRUE)

mn2beta(m, n, drop = TRUE)

# S3 method for class 'betaMix'
print(x, ...)

# S3 method for class 'betaBinomialMix'
print(x, ...)

# S3 method for class 'betaMix'
summary(object, probs = c(0.025, 0.5, 0.975), ...)

# S3 method for class 'betaBinomialMix'
summary(object, probs = c(0.025, 0.5, 0.975), ...)
```

## Arguments

- ...:

  List of mixture components.

- param:

  Determines how the parameters in the list are interpreted. See
  details.

- m:

  Vector of means of beta mixture components.

- s:

  Vector of standard deviations of beta mixture components.

- drop:

  Delete the dimensions of an array which have only one level.

- n:

  Vector of number of observations.

- x:

  The mixture to print

- object:

  Beta mixture object.

- probs:

  Quantiles reported by the `summary` function.

## Value

`mixbeta` returns a beta mixture with the specified mixture components.
`ms2beta` and `mn2beta` return the equivalent natural `a` and `b`
parametrization given parameters `m`, `s`, or `n`.

## Details

Each entry in the `...` argument list is expected to be a triplet of
numbers which defines the weight \\w_k\\, first and second parameter of
the mixture component \\k\\. A triplet can optionally be named which
will be used appropriately.

The first and second parameter can be given in different
parametrizations which is set by the `param` option:

- ab:

  Natural parametrization of Beta density (`a`=shape1 and `b`=shape2).
  Default.

- ms:

  Mean and standard deviation, \\m=a/(a+b)\\ and
  \\s=\sqrt{\frac{m(1-m)}{1+n}}\\, where \\n=a+b\\ is the number of
  observations. Note that \\s\\ must be less than \\\sqrt{m(1-m)}\\.

- mn:

  Mean and number of observations, \\n=a+b\\.

## See also

Other mixdist:
[`mix`](https://opensource.nibr.com/RBesT/reference/mix.md),
[`mixcombine()`](https://opensource.nibr.com/RBesT/reference/mixcombine.md),
[`mixgamma()`](https://opensource.nibr.com/RBesT/reference/mixgamma.md),
[`mixjson`](https://opensource.nibr.com/RBesT/reference/mixjson.md),
[`mixmvnorm()`](https://opensource.nibr.com/RBesT/reference/mixmvnorm.md),
[`mixnorm()`](https://opensource.nibr.com/RBesT/reference/mixnorm.md),
[`mixplot`](https://opensource.nibr.com/RBesT/reference/mixplot.md)

## Examples

``` r
## a beta mixture
bm <- mixbeta(rob = c(0.2, 2, 10), inf = c(0.4, 10, 100), inf2 = c(0.4, 30, 80))

# mean/standard deviation parametrization
bm2 <- mixbeta(rob = c(0.2, 0.3, 0.2), inf = c(0.8, 0.4, 0.01), param = "ms")

# mean/observations parametrization
bm3 <- mixbeta(rob = c(0.2, 0.3, 5), inf = c(0.8, 0.4, 30), param = "mn")

# even mixed is possible
bm4 <- mixbeta(rob = c(0.2, mn2beta(0.3, 5)), inf = c(0.8, ms2beta(0.4, 0.1)))

# print methods are defined
bm4
#> Univariate beta mixture
#> Mixture Components:
#>   rob  inf 
#> w  0.2  0.8
#> a  1.5  9.2
#> b  3.5 13.8
print(bm4)
#> Univariate beta mixture
#> Mixture Components:
#>   rob  inf 
#> w  0.2  0.8
#> a  1.5  9.2
#> b  3.5 13.8
```
