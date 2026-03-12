# Find root of univariate function of integers

Uses a bisectioning algorithm to search the give interval for a change
of sign and returns the integer which is closest to 0.

## Usage

``` r
uniroot_int(
  f,
  interval,
  ...,
  f.lower = f(interval[1], ...),
  f.upper = f(interval[2], ...),
  maxIter = 1000
)
```
