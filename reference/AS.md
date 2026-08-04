# Ankylosing Spondylitis.

Data set containing historical information for placebo for a phase II
trial of ankylosing spondylitis patients. The primary efficacy endpoint
was the percentage of patients with a 20% response according to the
Assessment of SpondyloArthritis international Society criteria for
improvement (ASAS20) at week 6.

## Usage

``` r
AS
```

## Format

A data frame with 8 rows and 3 variables:

- study:

  study

- n:

  study size

- r:

  number of events

## References

Baeten D, others (2013). “Anti-interleukin-17A monoclonal antibody
secukinumab in treatment of ankylosing spondylitis: a randomised,
double-blind, placebo-controlled trial.” *The Lancet*, **382**(9906),
1705–1713.

## Examples

``` r
.user_mc_options <- options()

set.seed(34563)
map_AS <- gMAP(cbind(r, n - r) ~ 1 | study,
  family = binomial,
  data = AS,
  tau.dist = "HalfNormal", tau.prior = 1,
  beta.prior = 2
)
#> Assuming default prior location   for beta: 0
## Recover user set sampling defaults
options(.user_mc_options)
```
