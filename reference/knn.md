# k nearest neighbor algorithm for multi-variate data

k nearest neighbor algorithm for multi-variate data

## Usage

``` r
knn(X, K = 2, init, Ninit = 50, verbose = FALSE, tol, Niter.max = 500)
```

## Arguments

- X:

  data matrix, i.e. observations X dimensions

- K:

  number of clusters to use

- init:

  list of p and mu used for initialization

- Ninit:

  number of samples used per cluster if no init argument is given

- verbose:

  allows print out of progress information; in verbose mode the cluster
  memberships are added to the output

- tol:

  smaller changes than tol in the objective function indicate
  convergence, if missing chosen automatically to be 1/5 of the smallest
  sample variance per dimension

- Niter.max:

  maximum number of admissible iterations
