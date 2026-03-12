# Decision Boundary for 1 Sample Designs

Calculates the decision boundary for a 1 sample design. This is the
critical value at which the decision function will change from 0
(failure) to 1 (success).

## Usage

``` r
decision1S_boundary(prior, n, decision, ...)

# S3 method for class 'betaMix'
decision1S_boundary(prior, n, decision, ...)

# S3 method for class 'normMix'
decision1S_boundary(prior, n, decision, sigma, eps = 1e-06, ...)

# S3 method for class 'gammaMix'
decision1S_boundary(prior, n, decision, eps = 1e-06, ...)
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

Returns the critical value \\y_c\\. For two-sided decision functions a
named vector with components `lower_or_equal_than` and `higher_than` is
returned, containing the critical values for the lower and upper
decision boundaries.

## Details

The specification of the 1 sample design (prior, sample size and
decision function, \\D(y)\\), uniquely defines the decision boundary

\$\$y_c = \max_y\\D(y) = 1\\,\$\$

which is the maximal value of \\y\\ whenever the decision \\D(y)\\
function changes its value from 1 to 0 for a decision function with
`lower.tail=TRUE` (otherwise the definition is \\y_c = \max\_{y}\\D(y) =
0\\\\). The decision function may change at most at a single critical
value as only one-sided decision functions are supported. Here, \\y\\ is
defined for binary and Poisson endpoints as the sufficient statistic \\y
= \sum\_{i=1}^{n} y_i\\ and for the normal case as the mean \\\bar{y} =
1/n \sum\_{i=1}^n y_i\\.

The convention for the critical value \\y_c\\ depends on whether a left
(`lower.tail=TRUE`) or right-sided decision function
(`lower.tail=FALSE`) is used. For `lower.tail=TRUE` the critical value
\\y_c\\ is the largest value for which the decision is 1, \\D(y \leq
y_c) = 1\\, while for `lower.tail=FALSE` then \\D(y \> y_c) = 1\\ holds.
This is aligned with the cumulative density function definition within R
(see for example [`pbinom()`](https://rdrr.io/r/stats/Binomial.html)).

## Methods (by class)

- `decision1S_boundary(betaMix)`: Applies for binomial model with a
  mixture beta prior. The calculations use exact expressions.

- `decision1S_boundary(normMix)`: Applies for the normal model with
  known standard deviation \\\sigma\\ and a normal mixture prior for the
  mean. As a consequence from the assumption of a known standard
  deviation, the calculation discards sampling uncertainty of the second
  moment. The function `decision1S_boundary` has an extra argument `eps`
  (defaults to \\10^{-6}\\). The critical value \\y_c\\ is searched in
  the region of probability mass `1-eps` for \\y\\.

- `decision1S_boundary(gammaMix)`: Applies for the Poisson model with a
  gamma mixture prior for the rate parameter. The function
  `decision1S_boundary` takes an extra argument `eps` (defaults to
  \\10^{-6}\\) which determines the region of probability mass `1-eps`
  where the boundary is searched for \\y\\.

## See also

Other design1S:
[`decision1S()`](https://opensource.nibr.com/RBesT/reference/decision1S.md),
[`oc1S()`](https://opensource.nibr.com/RBesT/reference/oc1S.md),
[`pos1S()`](https://opensource.nibr.com/RBesT/reference/pos1S.md)

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

# double criterion design
# statistical significance (like NI design)
dec1 <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)
# require mean to be at least as good as theta_c
dec2 <- decision1S(0.5, theta_c, lower.tail = TRUE)
# combination
decComb <- decision1S(c(1 - alpha, 0.5), c(theta_ni, theta_c), lower.tail = TRUE)

# critical value of double criterion design
decision1S_boundary(flat_prior, nL, decComb)
#> Using default prior reference scale 2
#> [1] 0.1357511

# ... is limited by the statistical significance ...
decision1S_boundary(flat_prior, nL, dec1)
#> Using default prior reference scale 2
#> [1] 0.1844494

# ... or the indecision point (whatever is smaller)
decision1S_boundary(flat_prior, nL, dec2)
#> Using default prior reference scale 2
#> [1] 0.1357511
```
