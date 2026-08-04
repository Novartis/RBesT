# Transform `gMAP` to `draws` objects

Transform a `gMAP` object to a format supported by the posterior
package.

## Usage

``` r
# S3 method for class 'gMAP'
as_draws(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'gMAP'
as_draws_matrix(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'gMAP'
as_draws_array(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'gMAP'
as_draws_df(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'gMAP'
as_draws_list(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'gMAP'
as_draws_rvars(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)
```

## Arguments

- x:

  A `gMAP` object.

- variable:

  A character vector providing the variables to extract. By default, all
  variables are extracted.

- regex:

  Logical; Should variable be treated as a (vector of) regular
  expressions? Any variable in `x` matching at least one of the regular
  expressions will be selected. Defaults to `FALSE`.

- inc_warmup:

  Should warmup draws be included? Defaults to `FALSE`.

- ...:

  Arguments passed to individual methods (if applicable).

## Details

To subset iterations, chains, or draws, use the
[`posterior::subset_draws()`](https://mc-stan.org/posterior/reference/subset_draws.html)
method after transforming the input object to a `draws` object.

## See also

[`posterior::draws()`](https://mc-stan.org/posterior/reference/draws.html)
[`posterior::subset_draws()`](https://mc-stan.org/posterior/reference/subset_draws.html)

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

post_AS <- as_draws(map_AS)

## Recover user set sampling defaults
options(.user_mc_options)
```
