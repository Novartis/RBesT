# Write and Read a Mixture Object with JSON

**\[experimental\]**

These functions write and read a mixture object in the JSON format.

## Usage

``` r
write_mix_json(mix, con, ...)

read_mix_json(con, ..., rescale = TRUE)
```

## Arguments

- mix:

  A mixture object to be saved to JSON.

- con:

  A connection specifying where the JSON will be written to or read.

- ...:

  Additional arguments passed to the
  [`jsonlite::toJSON()`](https://jeroen.r-universe.dev/jsonlite/reference/fromJSON.html)
  and
  [`jsonlite::fromJSON()`](https://jeroen.r-universe.dev/jsonlite/reference/fromJSON.html)
  function for writing and reading, respectively.

- rescale:

  A logical value indicating whether to rescale the mixture weights so
  that they sum to 1. Defaults to `TRUE`.

## Value

The `write_mix_json` function does not return a value while the
`read_mix_json` returns the mixture object stored in the connection
specified.

## Details

The mixture objects are written or read from the connection `con`, which
can be a character string specifying a file path or a connection object
as detailed in
[`base::connections()`](https://rdrr.io/r/base/connections.html).

When writing mixture objects as JSON it is strongly recommended to
explicitly set the number of digits (argument `digits`) to be used for
the numerical representation in order to control the accuracy of the
JSON representation of the mixture object. If the mixture object
inherits from the `"EM"` class (as is the case when the mixture is
created using the
[`mixfit()`](https://opensource.nibr.com/RBesT/reference/mixfit.md)
function), then the mixture object will be cast to a simple mixture
object such that diagnostics from the `"EM"` fitting procedure are
dropped from the object. For easier readability the user is encouraged
to set the argument `pretty=TRUE`, which is passed to the
[`jsonlite::toJSON()`](https://jeroen.r-universe.dev/jsonlite/reference/fromJSON.html)
function and makes the output more human readable.

Note that when reading in mixture objects, then these are not
necessarily equal to the mixtures passed to the `write_mix_json`
function. This is a consequence of the limited precision of the textual
representation as defined by the `digits` argument.

## See also

Other mixdist:
[`mix`](https://opensource.nibr.com/RBesT/reference/mix.md),
[`mixbeta()`](https://opensource.nibr.com/RBesT/reference/mixbeta.md),
[`mixcombine()`](https://opensource.nibr.com/RBesT/reference/mixcombine.md),
[`mixgamma()`](https://opensource.nibr.com/RBesT/reference/mixgamma.md),
[`mixmvnorm()`](https://opensource.nibr.com/RBesT/reference/mixmvnorm.md),
[`mixnorm()`](https://opensource.nibr.com/RBesT/reference/mixnorm.md),
[`mixplot`](https://opensource.nibr.com/RBesT/reference/mixplot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)

write_mix_json(nm, "normal_mixture.json", pretty=TRUE, digits=1)

mix <- read_mix_json("normal_mixture.json")
} # }
```
