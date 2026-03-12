# Probability of Success for 2 Sample Design

The `pos2S` function defines a 2 sample design (priors, sample sizes &
decision function) for the calculation of the probability of success. A
function is returned which calculates the calculates the frequency at
which the decision function is evaluated to 1 when parameters are
distributed according to the given distributions.

## Usage

``` r
pos2S(prior1, prior2, n1, n2, decision, ...)

# S3 method for class 'betaMix'
pos2S(prior1, prior2, n1, n2, decision, eps, ...)

# S3 method for class 'normMix'
pos2S(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps = 1e-06,
  Ngrid = 10,
  ...
)

# S3 method for class 'gammaMix'
pos2S(prior1, prior2, n1, n2, decision, eps = 1e-06, ...)
```

## Arguments

- prior1:

  Prior for sample 1.

- prior2:

  Prior for sample 2.

- n1, n2:

  Sample size of the respective samples. Sample size `n1` must be
  greater than 0 while sample size `n2` must be greater or equal to 0.

- decision:

  Two-sample decision function to use; see
  [`decision2S`](https://opensource.nibr.com/RBesT/reference/decision2S.md).

- ...:

  Optional arguments.

- eps:

  Support of random variables are determined as the interval covering
  `1-eps` probability mass. Defaults to \\10^{-6}\\.

- sigma1:

  The fixed reference scale of sample 1. If left unspecified, the
  default reference scale of the prior 1 is assumed.

- sigma2:

  The fixed reference scale of sample 2. If left unspecified, the
  default reference scale of the prior 2 is assumed.

- Ngrid:

  Determines density of discretization grid on which decision function
  is evaluated (see below for more details).

## Value

Returns a function which when called with two arguments `mix1` and
`mix2` will return the frequencies at which the decision function is
evaluated to 1. Each argument is expected to be a mixture distribution
representing the assumed true distribution of the parameter in each
group.

## Details

The `pos2S` function defines a 2 sample design and returns a function
which calculates its probability of success. The probability of success
is the frequency with which the decision function is evaluated to 1
under the assumption of a given true distribution of the data implied by
a distirbution of the parameters \\\theta_1\\ and \\\theta_2\\.

The calculation is analogous to the operating characeristics
[`oc2S()`](https://opensource.nibr.com/RBesT/reference/oc2S.md) with the
difference that instead of assuming known (point-wise) true parameter
values a distribution is specified for each parameter.

Calling the `pos2S` function calculates the decision boundary
\\D_1(y_2)\\ and returns a function which can be used to evaluate the
PoS for different predictive distributions. It is evaluated as

\$\$ \int\int\int f_2(y_2\|\theta_2) \\ p(\theta_2) \\
F_1(D_1(y_2)\|\theta_1) \\ p(\theta_1) \\ dy_2 d\theta_2 d\theta_1. \$\$

where \\F\\ is the distribution function of the sampling distribution
and \\p(\theta_1)\\ and \\p(\theta_2)\\ specifies the assumed true
distribution of the parameters \\\theta_1\\ and \\\theta_2\\,
respectively. Each distribution \\p(\theta_1)\\ and \\p(\theta_2)\\ is a
mixture distribution and given as the `mix1` and `mix2` argument to the
function.

For example, in the binary case an integration of the predictive
distribution, the BetaBinomial, instead of the binomial distribution
will be performed over the data space wherever the decision function is
evaluated to 1. All other aspects of the calculation are as for the
2-sample operating characteristics, see
[`oc2S()`](https://opensource.nibr.com/RBesT/reference/oc2S.md).

## Methods (by class)

- `pos2S(betaMix)`: Applies for binomial model with a mixture beta
  prior. The calculations use exact expressions. If the optional
  argument `eps` is defined, then an approximate method is used which
  limits the search for the decision boundary to the region of `1-eps`
  probability mass. This is useful for designs with large sample sizes
  where an exact approach is very costly to calculate.

- `pos2S(normMix)`: Applies for the normal model with known standard
  deviation \\\sigma\\ and normal mixture priors for the means. As a
  consequence from the assumption of a known standard deviation, the
  calculation discards sampling uncertainty of the second moment. The
  function has two extra arguments (with defaults): `eps` (\\10^{-6}\\)
  and `Ngrid` (10). The decision boundary is searched in the region of
  probability mass `1-eps`, respectively for \\y_1\\ and \\y_2\\. The
  continuous decision function is evaluated at a discrete grid, which is
  determined by a spacing with \\\delta_2 = \sigma_2/\sqrt{N\_{grid}}\\.
  Once the decision boundary is evaluated at the discrete steps, a
  spline is used to inter-polate the decision boundary at intermediate
  points.

- `pos2S(gammaMix)`: Applies for the Poisson model with a gamma mixture
  prior for the rate parameter. The function `pos2S` takes an extra
  argument `eps` (defaults to \\10^{-6}\\) which determines the region
  of probability mass `1-eps` where the boundary is searched for \\y_1\\
  and \\y_2\\, respectively.

## See also

Other design2S:
[`decision2S()`](https://opensource.nibr.com/RBesT/reference/decision2S.md),
[`decision2S_boundary()`](https://opensource.nibr.com/RBesT/reference/decision2S_boundary.md),
[`oc2S()`](https://opensource.nibr.com/RBesT/reference/oc2S.md)

## Examples

``` r
# see ?decision2S for details of example
priorT <- mixnorm(c(1, 0, 0.001), sigma = 88, param = "mn")
priorP <- mixnorm(c(1, -49, 20), sigma = 88, param = "mn")
# the success criteria is for delta which are larger than some
# threshold value which is why we set lower.tail=FALSE
successCrit <- decision2S(c(0.95, 0.5), c(0, 50), FALSE)

# example interim outcome
postP_interim <- postmix(priorP, n = 10, m = -50)
#> Using default prior reference scale 88
postT_interim <- postmix(priorT, n = 20, m = -80)
#> Using default prior reference scale 88

# assume that mean -50 / -80 were observed at the interim for
# placebo control(n=10) / active treatment(n=20) which gives
# the posteriors
postP_interim
#> Univariate normal mixture
#> Reference scale: 88
#> Mixture Components:
#>   comp1    
#> w   1.00000
#> m -49.33333
#> s  16.06653
postT_interim
#> Univariate normal mixture
#> Reference scale: 88
#> Mixture Components:
#>   comp1    
#> w   1.00000
#> m -79.99600
#> s  19.67691

# then the PoS to succeed after another 20/30 patients is
pos_final <- pos2S(postP_interim, postT_interim, 20, 30, successCrit)
#> Using default prior 1 reference scale 88
#> Using default prior 2 reference scale 88

pos_final(postP_interim, postT_interim)
#> [1] 0.145567
```
