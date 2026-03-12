# Exact Confidence interval for Binary Proportion

This function calculates the exact confidendence interval for a response
rate presented by \\n\\ and \\r\\.

## Usage

``` r
BinaryExactCI(r, n, alpha = 0.05, drop = TRUE)
```

## Arguments

- r:

  Number of success or responder

- n:

  Sample size

- alpha:

  confidence level

- drop:

  Determines if [`drop()`](https://rdrr.io/r/base/drop.html) will be
  called on the result

## Value

100 (1-\\\alpha\\)\\ response rate

## Details

Confidence intervals are obtained by a procedure first given in Clopper
and Pearson (1934). This guarantees that the confidence level is at
least (1-\\\alpha\\).

Details can be found in the publication listed below.

## References

Clopper, C. J. & Pearson, E. S. The use of confidence or fiducial limits
illustrated in the case of the binomial. Biometrika 1934.

## Examples

``` r
BinaryExactCI(3, 20, 0.05)
#>       2.5%      97.5% 
#> 0.03207094 0.37892683 
```
