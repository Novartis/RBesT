# Transform Densities with a link function

One-to-one transforms (mixture) of densities using a link function.

## Usage

``` r
dlink(object) <- value
```

## Arguments

- object:

  Mixture density to apply link to.

- value:

  Link.

  Note: link functions are assumed to be order preserving, i.e. if x_1
  \< x_2 holds, then link(x_1) \< link(x_2).
