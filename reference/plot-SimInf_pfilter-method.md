# Diagnostic plot of a particle filter object

Diagnostic plot of a particle filter object

## Usage

``` r
# S4 method for class 'SimInf_pfilter'
plot(x, y, ...)
```

## Arguments

- x:

  The `SimInf_pfilter` object to plot.

- y:

  If y is `NULL` or missing (default), the filtered trajectory (top) and
  the effective sample size (bottom) are displayed. If `y` is a
  character vector or a formula, the plot function for a `SimInf_model`
  object is called with the filtered trajectory, see
  [`plot,SimInf_model-method`](http://stewid.github.io/SimInf/reference/plot.md)
  for more details about the specification a plot.

- ...:

  Other graphical parameters that are passed on to the plot function.
