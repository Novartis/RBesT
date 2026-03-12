# Decision Function for 2 Sample Designs

The function sets up a 2 sample decision function with an arbitrary
number of conditions on the difference distribution.

## Usage

``` r
decision2S(
  pc = 0.975,
  qc = 0,
  lower.tail = TRUE,
  link = c("identity", "logit", "log")
)

oc2Sdecision(
  pc = 0.975,
  qc = 0,
  lower.tail = TRUE,
  link = c("identity", "logit", "log")
)
```

## Arguments

- pc:

  Vector of critical cumulative probabilities of the difference
  distribution.

- qc:

  Vector of respective critical values of the difference distribution.
  Must match the length of `pc`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(X \leq x)\\,
  otherwise, \\P(X \> x)\\.

- link:

  Enables application of a link function prior to evaluating the
  difference distribution. Can take one of the values `identity`
  (default), `logit` or `log`.

## Value

The function returns a decision function, of class `decision2S_1sided`
for one-sided, and of class `decision2S_2sided` for two-sided decisions.

One-sided decision functions take three arguments. The first and second
argument are expected to be mixture (posterior) distributions from which
the difference distribution is formed and all conditions are tested. The
third argument determines if the function acts as an indicator function
or if the function returns the distance from the decision boundary for
each condition in log-space. That is, the distance is 0 at the decision
boundary, negative for a 0 decision and positive for a 1 decision.

For two-sided decision functions, the two components can be extracted
with functions
[`lower()`](https://opensource.nibr.com/RBesT/reference/decision1S.md)
and
[`upper()`](https://opensource.nibr.com/RBesT/reference/decision1S.md).
The distance as calculated by the decision function is returned as a
list with components `lower` and `upper`.

## Details

This function creates a one- or two-sided decision function on the basis
of the difference distribution in a 2 sample situation. To support
double criterion designs, see *Neuenschwander et al., 2010*, an
arbitrary number of criterions can be given. The decision function
demands that the probability mass below the critical value `qc` of the
difference \\\theta_1 - \theta_2\\ is at least `pc`. Hence, for
`lower.tail=TRUE` condition \\i\\ is equivalent to

\$\$P(\theta_1 - \theta_2 \leq q\_{c,i}) \> p\_{c,i}\$\$

and the decision function is implemented as indicator function using the
heavy-side step function \\H(x)\\ which is \\0\\ for \\x \leq 0\\ and
\\1\\ for \\x \> 0\\. As all conditions must be met, the final indicator
function returns

\$\$\Pi_i H_i(P(\theta_1 - \theta_2 \leq q\_{c,i}) - p\_{c,i} ),\$\$

which is \\1\\ if all conditions are met and \\0\\ otherwise. For
`lower.tail=FALSE` differences must be greater than the given quantiles
`qc`.

For the case of a boolean vector given to `lower.tail` the direction of
each decision aligns respectively, and a two-sided decision function is
created.

Note that whenever a `link` other than `identity` is requested, then the
underlying densities are first transformed using the link function and
then the probabilties for the differences are calculated in the
transformed space. Hence, for a binary endpoint the default `identity`
link will calculate risk differences, the `logit` link will lead to
decisions based on the differences in `logit`s corresponding to a
criterion based on the log-odds. The `log` link will evaluate ratios
instead of absolute differences which could be useful for a binary
endpoint or counting rates. The respective critical quantiles `qc` must
be given on the transformed scale.

## Functions

- `oc2Sdecision()`: **\[deprecated\]** Deprecated old function name.
  Please use `decision2S` instead.

## References

Gsponer T, Gerber F, Bornkamp B, Ohlssen D, Vandemeulebroecke M,
Schmidli H.A practical guide to Bayesian group sequential designs.
*Pharm. Stat.*. 2014; 13: 71-80

## See also

Other design2S:
[`decision2S_boundary()`](https://opensource.nibr.com/RBesT/reference/decision2S_boundary.md),
[`oc2S()`](https://opensource.nibr.com/RBesT/reference/oc2S.md),
[`pos2S()`](https://opensource.nibr.com/RBesT/reference/pos2S.md)

## Examples

``` r
# see Gsponer et al., 2010
priorT <- mixnorm(c(1, 0, 0.001), sigma = 88, param = "mn")
priorP <- mixnorm(c(1, -49, 20), sigma = 88, param = "mn")
# the success criteria is for delta which are larger than some
# threshold value which is why we set lower.tail=FALSE
successCrit <- decision2S(c(0.95, 0.5), c(0, 50), FALSE)
# the futility criterion acts in the opposite direction
futilityCrit <- decision2S(c(0.90), c(40), TRUE)
# intermediate criteria can also be defined: no futility and no success.
# version 1: not significant, but good effect size.
intermediateCrit1 <- decision2S(
  c(1 - 0.95, 0.5, 1 - 0.90),
  c(0, 50, 40),
  c(TRUE, FALSE, TRUE)
)
# version 2: significant, but too small effect size.
intermediateCrit2 <- decision2S(
  c(0.95, 1 - 0.5, 1 - 0.90),
  c(0, 50, 40),
  c(FALSE, TRUE, TRUE)
)
# version 3: not significant and small effect size.
intermediateCrit3 <- decision2S(
  c(1 - 0.95, 1 - 0.5, 1 - 0.90),
  c(0, 50, 40),
  c(TRUE, TRUE, TRUE)
)

print(successCrit)
#> 2 sample decision function
#> Conditions for acceptance:
#> P(theta1 - theta2 > 0) > 0.95
#> P(theta1 - theta2 > 50) > 0.5
#> Link: identity 
print(futilityCrit)
#> 2 sample decision function
#> Conditions for acceptance:
#> P(theta1 - theta2 <= 40) > 0.9
#> Link: identity 
print(intermediateCrit1)
#> 2 sample two-sided decision function
#> Lower side conditions for acceptance:
#> P(theta1 - theta2 <= 0) > 0.05
#> P(theta1 - theta2 <= 40) > 0.1
#> Upper side conditions for acceptance:
#> P(theta1 - theta2 > 50) > 0.5
#> Link: identity 
print(intermediateCrit2)
#> 2 sample two-sided decision function
#> Lower side conditions for acceptance:
#> P(theta1 - theta2 <= 50) > 0.5
#> P(theta1 - theta2 <= 40) > 0.1
#> Upper side conditions for acceptance:
#> P(theta1 - theta2 > 0) > 0.95
#> Link: identity 
print(intermediateCrit3)
#> 2 sample decision function
#> Conditions for acceptance:
#> P(theta1 - theta2 <= 0) > 0.05
#> P(theta1 - theta2 <= 50) > 0.5
#> P(theta1 - theta2 <= 40) > 0.1
#> Link: identity 

# consider decision for specific outcomes
postP_interim <- postmix(priorP, n = 10, m = -50)
#> Using default prior reference scale 88
postT_interim <- postmix(priorT, n = 20, m = -80)
#> Using default prior reference scale 88
futilityCrit(postP_interim, postT_interim)
#> [1] 0
successCrit(postP_interim, postT_interim)
#> [1] 0
intermediateCrit1(postP_interim, postT_interim)
#> [1] 0
intermediateCrit2(postP_interim, postT_interim)
#> [1] 0
intermediateCrit3(postP_interim, postT_interim)
#> [1] 1

# Binary endpoint with double criterion decision on log-odds scale
# 95% certain positive difference and an odds ratio of 2 at least
decL2 <- decision2S(c(0.95, 0.5), c(0, log(2)), lower.tail = FALSE, link = "logit")
# 95% certain positive difference and an odds ratio of 3 at least
decL3 <- decision2S(c(0.95, 0.5), c(0, log(3)), lower.tail = FALSE, link = "logit")

# data scenario
post1 <- postmix(mixbeta(c(1, 1, 1)), n = 40, r = 10)
post2 <- postmix(mixbeta(c(1, 1, 1)), n = 40, r = 18)

# positive outcome and a median odds ratio of at least 2 ...
decL2(post2, post1)
#> [1] 1
# ... but not more than 3
decL3(post2, post1)
#> [1] 0
```
