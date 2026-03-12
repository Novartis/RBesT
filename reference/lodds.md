# Logit (log-odds) and inverse-logit function.

Calculates the logit (log-odds) and inverse-logit.

## Usage

``` r
logit(mu)

inv_logit(eta)
```

## Arguments

- mu:

  A numeric object with probabilies, with values in the in the range
  \\\[0,1\]\\. Missing values (`NA`s) are allowed.

- eta:

  A numeric object with log-odds values, with values in the range
  \\\[-\infty,\infty\]\\. Missing values (`NA`s) are allowed.

## Value

A numeric object of the same type as mu and eta containing the logits or
inverse logit of the input values. The logit and inverse transformation
equates to

\$\$\text{logit}(\mu) = \log(\mu/(1-\mu))\$\$
\$\$\text{logit}^{-1}(\eta)= \exp(\eta)/(1 + \exp(\eta)).\$\$

## Details

Values of mu equal to 0 or 1 will return \\-\infty\\ or \\\infty\\
respectively.

## Examples

``` r
logit(0.2)
#> [1] -1.386294
inv_logit(-1.386)
#> [1] 0.2000471
```
