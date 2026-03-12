# Probability of Success for a 1 Sample Design

The `pos1S` function defines a 1 sample design (prior, sample size,
decision function) for the calculation of the frequency at which the
decision is evaluated to 1 when assuming a distribution for the
parameter. A function is returned which performs the actual operating
characteristics calculations.

## Usage

``` r
pos1S(prior, n, decision, ...)

# S3 method for class 'betaMix'
pos1S(prior, n, decision, ...)

# S3 method for class 'normMix'
pos1S(prior, n, decision, sigma, eps = 1e-06, ...)

# S3 method for class 'gammaMix'
pos1S(prior, n, decision, eps = 1e-06, ...)
```

## Arguments

- prior:

  Prior for analysis.

- n:

  Sample size for the experiment.

- decision:

  One-sample decision function to use; see
  [`decision1S`](https://opensource.nibr.com/RBesT/reference/decision1S.md).

- ...:

  Optional arguments.

- sigma:

  The fixed reference scale. If left unspecified, the default reference
  scale of the prior is assumed.

- eps:

  Support of random variables are determined as the interval covering
  `1-eps` probability mass. Defaults to \\10^{-6}\\.

## Value

Returns a function that takes as single argument `mix`, which is the
mixture distribution of the control parameter. Calling this function
with a mixture distribution then calculates the PoS.

## Details

The `pos1S` function defines a 1 sample design and returns a function
which calculates its probability of success. The probability of success
is the frequency with which the decision function is evaluated to 1
under the assumption of a given true distribution of the data implied by
a distribution of the parameter \\\theta\\.

Calling the `pos1S` function calculates the critical value \\y_c\\ and
returns a function which can be used to evaluate the PoS for different
predictive distributions and is evaluated as

\$\$ \int F(y_c\|\theta) p(\theta) d\theta, \$\$

where \\F\\ is the distribution function of the sampling distribution
and \\p(\theta)\\ specifies the assumed true distribution of the
parameter \\\theta\\. The distribution \\p(\theta)\\ is a mixture
distribution and given as the `mix` argument to the function.

## Methods (by class)

- `pos1S(betaMix)`: Applies for binomial model with a mixture beta
  prior. The calculations use exact expressions.

- `pos1S(normMix)`: Applies for the normal model with known standard
  deviation \\\sigma\\ and a normal mixture prior for the mean. As a
  consequence from the assumption of a known standard deviation, the
  calculation discards sampling uncertainty of the second moment. The
  function `pos1S` has an extra argument `eps` (defaults to
  \\10^{-6}\\). The critical value \\y_c\\ is searched in the region of
  probability mass `1-eps` for \\y\\.

- `pos1S(gammaMix)`: Applies for the Poisson model with a gamma mixture
  prior for the rate parameter. The function `pos1S` takes an extra
  argument `eps` (defaults to \\10^{-6}\\) which determines the region
  of probability mass `1-eps` where the boundary is searched for \\y\\.

## See also

Other design1S:
[`decision1S()`](https://opensource.nibr.com/RBesT/reference/decision1S.md),
[`decision1S_boundary()`](https://opensource.nibr.com/RBesT/reference/decision1S_boundary.md),
[`oc1S()`](https://opensource.nibr.com/RBesT/reference/oc1S.md)

## Examples

``` r
# non-inferiority example using normal approximation of log-hazard
# ratio, see ?decision1S for all details
s <- 2
flat_prior <- mixnorm(c(1, 0, 100), sigma = s)
nL <- 233
theta_ni <- 0.4
theta_a <- 0
alpha <- 0.05
beta <- 0.2
za <- qnorm(1 - alpha)
zb <- qnorm(1 - beta)
n1 <- round((s * (za + zb) / (theta_ni - theta_a))^2)
theta_c <- theta_ni - za * s / sqrt(n1)

# assume we would like to conduct at an interim analysis
# of PoS after having observed 20 events with a HR of 0.8.
# We first need the posterior at the interim ...
post_ia <- postmix(flat_prior, m = log(0.8), n = 20)
#> Using default prior reference scale 2

# dual criterion
decComb <- decision1S(c(1 - alpha, 0.5), c(theta_ni, theta_c), lower.tail = TRUE)

# ... and we would like to know the PoS for a successful
# trial at the end when observing 10 more events
pos_ia <- pos1S(post_ia, 10, decComb)
#> Using default prior reference scale 2

# our knowledge at the interim is just the posterior at
# interim such that the PoS is
pos_ia(post_ia)
#> [1] 0.534741
```
