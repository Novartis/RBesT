# Diagnostic plots for gMAP analyses

Diagnostic plots for gMAP analyses

## Usage

``` r
# S3 method for class 'gMAP'
plot(x, size = NULL, linewidth = NULL, ...)
```

## Arguments

- x:

  [`gMAP()`](https://opensource.nibr.com/RBesT/reference/gMAP.md) object

- size:

  Controls size of forest plot.

- linewidth:

  Controls line sizes of traceplots.

- ...:

  Ignored.

## Value

The function returns a list of
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
objects.

## Details

Creates MCMC diagnostics and a forest plot (including model estimates)
for a [`gMAP()`](https://opensource.nibr.com/RBesT/reference/gMAP.md)
analysis. For a customized forest plot, please use the dedicated
function
[`forest_plot()`](https://opensource.nibr.com/RBesT/reference/forest_plot.md).

## Customizing ggplot2 plots

The returned plot is a ggplot2 object. Please refer to the "Customizing
Plots" vignette which is part of RBesT documentation for an
introduction. For simple modifications (change labels, add reference
lines, ...) consider the commands found in
[`bayesplot-helpers`](https://mc-stan.org/bayesplot/reference/bayesplot-helpers.html).
For more advanced customizations please use the ggplot2 package
directly. A description of the most common tasks can be found in the [R
Cookbook](http://www.cookbook-r.com/Graphs/) and a full reference of
available commands can be found at the [ggplot2 documentation
site](https://ggplot2.tidyverse.org/reference/).
