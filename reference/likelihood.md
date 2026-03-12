# Read and Set Likelihood to the Corresponding Conjugate Prior

Read and set the likelihood distribution corresponding to the conjugate
prior distribution.

## Usage

``` r
likelihood(mix)

likelihood(mix) <- value
```

## Arguments

- mix:

  Prior mixture distribution.

- value:

  New likelihood. **Should** only be changed for Gamma priors as these
  are supported with either Poisson (`value="poisson"`) or Exponential
  (`value="exp"`) likelihoods.

## Details

If the prior and posterior distributions are in the same family, then
the prior distribution is called a conjugate prior for the likelihood
function.

## Supported Conjugate Prior-Likelihood Pairs

|                     |                             |                             |                |
|---------------------|-----------------------------|-----------------------------|----------------|
| **Prior/Posterior** | **Likelihood**              | **Predictive**              | **Summaries**  |
| Beta                | Binomial                    | Beta-Binomial               | `n`, `r`       |
| Normal              | Normal (*fixed \\\sigma\\*) | Normal                      | `n`, `m`, `se` |
| Gamma               | Poisson                     | Gamma-Poisson               | `n`, `m`       |
| Gamma               | Exponential                 | Gamma-Exp (*not supported*) | `n`, `m`       |

## Examples

``` r
# Gamma mixture
gmix <- mixgamma(c(0.3, 20, 4), c(0.7, 50, 10))

# read out conjugate partner
likelihood(gmix)
#> [1] "poisson"

ess(gmix)
#> [1] 7.159378

# set conjugate partner
likelihood(gmix) <- "exp"

# ... which changes the interpretation of the mixture
ess(gmix)
#> [1] 34.93388
```
