# Resolve family/sigma arguments into a sigma_fun closure

Validates the family/sigma combination and returns a list with:

- `sigma_fun`: NULL (fixed-sigma path) or a closure `function(eta)`
  returning sigma at each eta.

- `family`: the family object (possibly set to NULL for
  gaussian(identity) short-circuit).

## Usage

``` r
resolve_sigma_family(family, sigma_missing, sigma = NULL, offset = 0)
```

## Arguments

- family:

  A family object or NULL.

- sigma_missing:

  Logical: was sigma missing in the calling function?

- sigma:

  Numeric sigma value (only used when sigma_missing is FALSE).

- offset:

  Numeric scalar offset.

## Value

A list with components `sigma_fun` and `family`.
