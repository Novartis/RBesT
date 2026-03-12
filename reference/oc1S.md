# Operating Characteristics for 1 Sample Design

The `oc1S` function defines a 1 sample design (prior, sample size,
decision function) for the calculation of the frequency at which the
decision is evaluated to 1 conditional on assuming known parameters. A
function is returned which performs the actual operating characteristics
calculations.

## Usage

``` r
oc1S(prior, n, decision, ...)

# S3 method for class 'betaMix'
oc1S(prior, n, decision, ...)

# S3 method for class 'normMix'
oc1S(prior, n, decision, sigma, eps = 1e-06, ...)

# S3 method for class 'gammaMix'
oc1S(prior, n, decision, eps = 1e-06, ...)
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

Returns a function with one argument `theta` which calculates the
frequency at which the decision function is evaluated to 1 for the
defined 1 sample design. Note that the returned function takes vectors
arguments.

## Details

The `oc1S` function defines a 1 sample design and returns a function
which calculates its operating characteristics. This is the frequency
with which the decision function is evaluated to 1 under the assumption
of a given true distribution of the data defined by a known parameter
\\\theta\\. The 1 sample design is defined by the prior, the sample size
and the decision function, \\D(y)\\. These uniquely define the decision
boundary, see
[`decision1S_boundary()`](https://opensource.nibr.com/RBesT/reference/decision1S_boundary.md).

When calling the `oc1S` function, then internally the critical value
\\y_c\\ (using
[`decision1S_boundary()`](https://opensource.nibr.com/RBesT/reference/decision1S_boundary.md))
is calculated and a function is returns which can be used to calculated
the desired frequency which is evaluated as

\$\$ F(y_c\|\theta). \$\$

## Methods (by class)

- `oc1S(betaMix)`: Applies for binomial model with a mixture beta prior.
  The calculations use exact expressions.

- `oc1S(normMix)`: Applies for the normal model with known standard
  deviation \\\sigma\\ and a normal mixture prior for the mean. As a
  consequence from the assumption of a known standard deviation, the
  calculation discards sampling uncertainty of the second moment. The
  function `oc1S` has an extra argument `eps` (defaults to \\10^{-6}\\).
  The critical value \\y_c\\ is searched in the region of probability
  mass `1-eps` for \\y\\.

- `oc1S(gammaMix)`: Applies for the Poisson model with a gamma mixture
  prior for the rate parameter. The function `oc1S` takes an extra
  argument `eps` (defaults to \\10^{-6}\\) which determines the region
  of probability mass `1-eps` where the boundary is searched for \\y\\.

## See also

Other design1S:
[`decision1S()`](https://opensource.nibr.com/RBesT/reference/decision1S.md),
[`decision1S_boundary()`](https://opensource.nibr.com/RBesT/reference/decision1S_boundary.md),
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

# standard NI design
decA <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)

# double criterion design
# statistical significance (like NI design)
dec1 <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)
# require mean to be at least as good as theta_c
dec2 <- decision1S(0.5, theta_c, lower.tail = TRUE)
# combination
decComb <- decision1S(c(1 - alpha, 0.5), c(theta_ni, theta_c), lower.tail = TRUE)

theta_eval <- c(theta_a, theta_c, theta_ni)

# evaluate different designs at two sample sizes
designA_n1 <- oc1S(flat_prior, n1, decA)
#> Using default prior reference scale 2
designA_nL <- oc1S(flat_prior, nL, decA)
#> Using default prior reference scale 2
designC_n1 <- oc1S(flat_prior, n1, decComb)
#> Using default prior reference scale 2
designC_nL <- oc1S(flat_prior, nL, decComb)
#> Using default prior reference scale 2

# evaluate designs at the key log-HR of positive treatment (HR<1),
# the indecision point and the NI margin

designA_n1(theta_eval)
#> [1] 0.80084529 0.49980774 0.04995032
designA_nL(theta_eval)
#> [1] 0.92039728 0.64489439 0.04997268
designC_n1(theta_eval)
#> [1] 0.80084529 0.49980774 0.04995032
designC_nL(theta_eval)
#> [1] 0.84991646 0.49995959 0.02185859

# to understand further the dual criterion design it is useful to
# evaluate the criterions separately:
# statistical significance criterion to warrant NI...
designC1_nL <- oc1S(flat_prior, nL, dec1)
#> Using default prior reference scale 2
# ... or the clinically determined indifference point
designC2_nL <- oc1S(flat_prior, nL, dec2)
#> Using default prior reference scale 2

designC1_nL(theta_eval)
#> [1] 0.92039728 0.64489439 0.04997268
designC2_nL(theta_eval)
#> [1] 0.84991646 0.49995959 0.02185859

# see also ?decision1S_boundary to see which of the two criterions
# will drive the decision

# it can also be of interest to evaluate intermediate decisions
# where the trial is significant for non-inferiority but
# the mean estimate is in an intermediate range, say between theta_c
# and theta_f:
theta_f <- 0.3
decCombIntermediate <- decision1S(
  c(1 - alpha, 0.5, 0.8),
  c(theta_ni, theta_c, theta_f),
  lower.tail = c(TRUE, FALSE, TRUE)
)

# evaluate at two sample sizes again:
designIntermediate_n1 <- oc1S(flat_prior, n1, decCombIntermediate)
#> Using default prior reference scale 2
designIntermediate_nL <- oc1S(flat_prior, nL, decCombIntermediate)
#> Using default prior reference scale 2

designIntermediate_n1(theta_eval)
#> [1] 0 0 0
designIntermediate_nL(theta_eval)
#> [1] 0.07043191 0.14485113 0.02810313
```
