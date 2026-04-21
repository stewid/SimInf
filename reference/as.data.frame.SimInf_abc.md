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

- **Parameter columns**: One column for each parameter estimated in the
  ABC analysis (e.g., `beta`, `gamma`, `sigma`). The column names match
  the parameter names defined in the `priors`.

## See also

[`abc`](http://stewid.github.io/SimInf/reference/abc.md) for running the
ABC analysis,
[`SimInf_abc`](http://stewid.github.io/SimInf/reference/SimInf_abc-class.md)
for the class definition, and
[`continue_abc`](http://stewid.github.io/SimInf/reference/continue_abc.md)
for continuing an existing ABC run.
