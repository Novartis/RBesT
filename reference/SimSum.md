# Summarize Arrays

The function calculates summary statistics from arbitrary arrays.

## Usage

``` r
SimSum(
  x,
  min.max = FALSE,
  n.sim = FALSE,
  probs = c(0.025, 0.5, 0.975),
  margin = ifelse(is.null(dim(x) | length(dim(x)) == 1), 2, length(dim(x)))
)
```

## Arguments

- x:

  Object to summarize which can be a numerical vector, matrix or a
  multi-dimensional array

- min.max:

  Enables to include minimum and maximum in the output.

- n.sim:

  Enables to include the number of observations in the output.

- probs:

  Quantiles to output.

- margin:

  Margin of the input array over which the summary function is applied.

## Details

The function calculates by default the mean, standard deviation and the
specified qantiles which are by default the median and the 95% interval.

If a mulit-dimensional array is specified as `x`, then the function will
by default calculate the summaries over the margin of the largest
dimension. For the case of a vector and a matrix, the function will
transpose the results for better readabiliy.
