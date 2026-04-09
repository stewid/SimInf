# Display the ABC posterior distribution

Produce diagnostic plots of the Approximate Bayesian Computation (ABC)
posterior distribution stored in a `SimInf_abc` object.

## Usage

``` r
# S4 method for class 'SimInf_abc'
plot(x, y, ...)
```

## Arguments

- x:

  The `SimInf_abc` object containing the ABC results.

- y:

  The generation number to plot. The default is `NULL`, which displays
  the **last** generation (the final posterior). Specify an integer to
  view intermediate generations for convergence diagnostics.

- ...:

  Additional graphical arguments passed to the underlying plotting
  functions (e.g., `col` for contour colors, `lwd`).

## Details

The function generates a scatterplot matrix of the parameter values for
the specified generation:

- **Diagonal panels**: Display a normalized density estimate with a rug
  plot showing individual samples.

- **Upper triangular panels**: Display scatterplots of the raw parameter
  samples for each pair of variables.

- **Lower triangular panels**: Display contour lines representing the 2D
  kernel density estimate of the joint distribution between parameter
  pairs.

If only a single parameter is selected, a single density plot with a rug
is produced.
