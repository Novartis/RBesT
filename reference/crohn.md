# Crohn's disease.

Data set containing historical information for placebo arm of relevant
studies for the treatment of Crohn's disease. The primary outcome is
change from baseline in Crohn's Disease Activity Index (CDAI) over a
duration of 6 weeks. Standard deviation of change from baseline endpoint
is approximately 88.

## Usage

``` r
crohn
```

## Format

A data frame with 4 rows and 3 variables:

- study:

  study

- n:

  study size

- y:

  mean CDAI change

## References

Hueber W, Sands BE, Lewitzky S, Vandemeulebroecke M, others (2012).
“Secukinumab, a human anti-IL-17A monoclonal antibody, for moderate to
severe Crohn's disease.” *Gut*, **61**(12), 1693–1700.

## Examples

``` r
.user_mc_options <- options()

set.seed(546346)
map_crohn <- gMAP(cbind(y, y.se) ~ 1 | study,
  family = gaussian,
  data = transform(crohn, y.se = 88 / sqrt(n)),
  weights = n,
  tau.dist = "HalfNormal", tau.prior = 44,
  beta.prior = cbind(0, 88)
)
## Recover user set sampling defaults
options(.user_mc_options)
```
