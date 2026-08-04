# Find root of univariate function of integers

Uses a bisectioning algorithm to search the give interval for a change
of sign and returns the integer which is closest to 0. When `extendInt`
is `"upX"` (f increasing through the root) or `"downX"` (f decreasing),
a bracket that does not enclose a sign change is widened outward
(doubling), bounded by `clamp`, until it brackets a root. If no sign
change exists within `clamp`, the constant sign is reported as `-Inf` (f
\< 0) or `Inf` (f \> 0); with `extendInt = "no"` the legacy empty
[`numeric()`](https://rdrr.io/r/base/numeric.html) is returned instead.

## Usage

``` r
uniroot_int(
  f,
  interval,
  ...,
  f.lower = f(interval[1], ...),
  f.upper = f(interval[2], ...),
  extendInt = c("no", "upX", "downX"),
  clamp = interval,
  maxIter = 1000
)
```
