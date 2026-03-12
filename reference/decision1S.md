# Decision Function for 1 Sample Designs

The function sets up a 1 sample decision function with an arbitrary
number of conditions.

## Usage

``` r
decision1S(pc = 0.975, qc = 0, lower.tail = TRUE)

has_lower(x)

has_upper(x)

lower(x)

upper(x)

oc1Sdecision(pc = 0.975, qc = 0, lower.tail = TRUE)
```

## Arguments

- pc:

  Vector of critical cumulative probabilities.

- qc:

  Vector of respective critical values. Must match the length of `pc`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(X \leq x)\\,
  otherwise, \\P(X \> x)\\. Either length 1 or same length as `pc`.

- x:

  Two-sided decision function.

## Value

The function returns a decision function (of class `decision1S_1sided`
for one-sided, and of class `decision1S_2sided` for two-sided decisions)
which takes two arguments. The first argument is expected to be a
mixture (posterior) distribution which is tested if the specified
conditions are met. The logical second argument determines if the
function acts as an indicator function or if the function returns the
distance from the decision boundary for each condition in log-space,
i.e. the distance is 0 at the decision boundary, negative for a 0
decision and positive for a 1 decision.

For two-sided decision functions, the two components can be extracted
with functions `lower()` and `upper()`. The distance as calculated by
the decision function is returned as a list with components `lower` and
`upper`.

## Details

For `lower.tail` being either `TRUE` or `FALSE`, the function creates a
one-sided decision function which takes two arguments. The first
argument is expected to be a mixture (posterior) distribution. This
distribution is tested whether it fulfills all the required threshold
conditions specified with the `pc` and `qc` arguments and returns 1 if
all conditions are met and 0 otherwise. Hence, for `lower.tail=TRUE`
condition \\i\\ is equivalent to

\$\$P(\theta \leq q\_{c,i}) \> p\_{c,i}\$\$

and the decision function is implemented as indicator function on the
basis of the heavy-side step function \\H(x)\\ which is \\0\\ for \\x
\leq 0\\ and \\1\\ for \\x \> 0\\. As all conditions must be met, the
final indicator function returns

\$\$\Pi_i H_i(P(\theta \leq q\_{c,i}) - p\_{c,i} ).\$\$

For the case of a boolean vector given to `lower.tail` the direction of
each decision aligns respectively, and a two-sided decision function is
created.

When the second argument is set to `TRUE` a distance metric is returned
component-wise per defined condition as

\$\$ D_i = \log(P(\theta \< q\_{c,i})) - \log(p\_{c,i}) .\$\$

These indicator functions can be used as input for 1-sample boundary, OC
or PoS calculations using
[`oc1S()`](https://opensource.nibr.com/RBesT/reference/oc1S.md) or
[`pos1S()`](https://opensource.nibr.com/RBesT/reference/pos1S.md) .

## Functions

- `oc1Sdecision()`: **\[deprecated\]** Deprecated old function name.
  Please use `decision1S` instead.

## References

Neuenschwander B, Rouyrre N, Hollaender H, Zuber E, Branson M. A proof
of concept phase II non-inferiority criterion. *Stat. in Med.*. 2011,
30:1618-1627

## See also

Other design1S:
[`decision1S_boundary()`](https://opensource.nibr.com/RBesT/reference/decision1S_boundary.md),
[`oc1S()`](https://opensource.nibr.com/RBesT/reference/oc1S.md),
[`pos1S()`](https://opensource.nibr.com/RBesT/reference/pos1S.md)

## Examples

``` r
# see Neuenschwander et al., 2011

# example is for a time-to-event trial evaluating non-inferiority (NI)
# using a normal approximation for the log-hazard ratio

# reference scale
s <- 2
theta_ni <- 0.4
theta_a <- 0
alpha <- 0.05
beta <- 0.2
za <- qnorm(1 - alpha)
zb <- qnorm(1 - beta)
n1 <- round((s * (za + zb) / (theta_ni - theta_a))^2) # n for which design was intended
nL <- 233
c1 <- theta_ni - za * s / sqrt(n1)

# flat prior
flat_prior <- mixnorm(c(1, 0, 100), sigma = s)

# standard NI design
decA <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)

# for double criterion with indecision point (mean estimate must be
# lower than this)
theta_c <- c1

# double criterion design
# statistical significance (like NI design)
dec1 <- decision1S(1 - alpha, theta_ni, lower.tail = TRUE)
# require mean to be at least as good as theta_c
dec2 <- decision1S(0.5, theta_c, lower.tail = TRUE)
# combination
decComb <- decision1S(c(1 - alpha, 0.5), c(theta_ni, theta_c), lower.tail = TRUE)

theta_eval <- c(theta_a, theta_c, theta_ni)

# we can display the decision function definition
decComb
#> 1 sample decision function
#> Conditions for acceptance:
#> P(theta <= 0.4) > 0.95
#> P(theta <= 0.13576435472344) > 0.5

# and use it to decide if a given distribution fulfills all
# criterions defined
# for the prior
decComb(flat_prior)
#> [1] 0
# or for a possible outcome of the trial
# here with HR of 0.8 for 40 events
decComb(postmix(flat_prior, m = log(0.8), n = 40))
#> Using default prior reference scale 2
#> [1] 1

# A two-sided decision function can be useful to determine if
# certain intermediate (i.e. neither "go" nor "stop") decisions
# are to be made based on the posterior distribution.
# For example, in the above situation we might have an intermediate
# scenario where the trial is significant for non-inferiority but
# the mean estimate is in an intermediate range, say between theta_c
# theta_f:
theta_f <- 0.3
decCombIntermediate <- decision1S(
  c(1 - alpha, 0.5, 0.8),
  c(theta_ni, theta_c, theta_f),
  lower.tail = c(TRUE, FALSE, TRUE)
)
# Not fulfilled for the prior:
decCombIntermediate(flat_prior)
#> [1] 0
# But for a hypothetical trial outcome with HR 1.2 and 300 events:
decCombIntermediate(postmix(flat_prior, m = log(1.2), n = 300))
#> Using default prior reference scale 2
#> [1] 1
```
