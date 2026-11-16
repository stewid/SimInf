# Display the ABC posterior distribution

Produce diagnostic plots of the Approximate Bayesian Computation (ABC)
posterior distribution stored in a `SimInf_abc` object. The parameter
values of the accepted particles in the selected generation are
displayed: with more than one parameter, a scatterplot matrix of the
pairwise distributions is produced, with one row and column per
parameter.

- **Diagonal panels**: Display a normalized density estimate with a rug
  plot showing individual samples. The density is estimated with the
  `bandwidth` `'SJ-ste'` and scaled so that its peak equals one.

- **Upper triangular panels**: Display scatterplots of the raw parameter
  samples for each pair of variables.

- **Lower triangular panels**: Display contour lines of the 2D kernel
  density estimate of the joint distribution between parameter pairs,
  computed with [`MASS::kde2d`](https://rdrr.io/pkg/MASS/man/kde2d.html)
  with default bandwidths and without contour labels.

If only a single parameter is selected, a single density plot with a rug
is produced.

## Usage

``` r
# S4 method for class 'SimInf_abc'
plot(x, y, ...)
```

## Arguments

- x:

  The `SimInf_abc` object containing the ABC results.

- y:

  The generation to plot as an integer in `1:n_generations(x)`. The
  default is `NULL`, which displays the last generation. Specify an
  earlier generation to inspect intermediate populations for convergence
  diagnostics.

- ...:

  Additional graphical arguments passed on to the underlying plotting
  functions, and therefore to the density estimates, contours, and
  scatterplots alike (e.g. `col`, `lwd`).
