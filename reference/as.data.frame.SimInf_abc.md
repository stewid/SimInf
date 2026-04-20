# Coerce a `SimInf_abc` object to a `data.frame`

Convert the results of an Approximate Bayesian Computation (ABC)
analysis into a single `data.frame`. This function extracts the particle
parameters, their acceptance weights, and the generation number for
every particle across all generations.

## Usage

``` r
# S3 method for class 'SimInf_abc'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `SimInf_abc` object.

- ...:

  Additional arguments (currently ignored).

## Value

A `data.frame` containing all particles from all generations.

## Details

The resulting `data.frame` has one row per particle. The columns
include:

- `generation`: The generation number (integer).

- `weight`: The normalized weight of the particle (numeric).

- `...`: Columns corresponding to the parameter names defined in the ABC
  analysis (e.g., `beta`, `gamma`).
