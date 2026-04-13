# Set a sparse recording template for simulation results

Define a template (or "punchcard") specifying which data points (nodes,
time-points, and compartments) should be recorded during a simulation.
This feature is useful for saving memory when the model has many nodes
or time-points, but only a subset of the data is needed for
post-processing.

## Usage

``` r
punchcard(model) <- value

# S4 method for class 'SimInf_model'
punchcard(model) <- value
```

## Arguments

- model:

  A `SimInf_model` object.

- value:

  A `data.frame` defining the recording template, or `NULL` to reset the
  model to record all compartments for all nodes at all time-points in
  `tspan`.

## Details

The template is specified as a `data.frame` with columns:

- `time`: The time-points to record.

- `node`: The node indices to record.

- `CompartmentName`: Logical columns (e.g., `S`, `I`, `R`) indicating
  whether to record that specific compartment for the corresponding row.
  If a compartment column is omitted, it defaults to `FALSE` (not
  recorded). If only `time` and `node` are provided, all compartments
  are recorded for those points.

**Important:** The specified `time` values, `node` indices, and
compartment names must exist in the model. If any value in the template
does not match the model's configuration (e.g., a time point not in
`tspan` or a compartment not defined in the model), an **error will be
raised** when the template is set.

Note that the template only affects which data is stored in the result
object; it does not affect how the solver simulates the trajectory.

## Examples

``` r
## For reproducibility, set the seed and number of threads.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 6 nodes
u0 <- data.frame(
  S = 100:105,
  I = 1:6,
  R = rep(0, 6)
)
model <- SIR(
  u0 = u0,
  tspan = 1:10,
  beta = 0.16,
  gamma = 0.077
)

## Run the model with default recording (all data)
result_full <- run(model)
head(trajectory(result_full))
#>   node time   S I R
#> 1    1    1 100 1 0
#> 2    2    1 101 2 0
#> 3    3    1 102 3 0
#> 4    4    1 102 5 0
#> 5    5    1 103 6 0
#> 6    6    1 105 6 0

## Define a template to record only nodes 2 and 4 at times 3 and 5
## for all compartments.
df <- data.frame(
  time = c(3, 5, 3, 5),
  node = c(2, 2, 4, 4),
  S = TRUE, I = TRUE, R = TRUE
)
punchcard(model) <- df

result_sparse <- run(model)
trajectory(result_sparse)
#>   node time   S I R
#> 1    2    3 100 3 0
#> 2    4    3 102 5 0
#> 3    2    5 100 3 0
#> 4    4    5 100 6 1

## Record only specific compartments (e.g., S and R, but not I)
df <- data.frame(
  time = c(3, 5, 3, 5),
  node = c(2, 2, 4, 4),
  S = TRUE, I = FALSE, R = TRUE
)
punchcard(model) <- df
result_partial <- run(model)
trajectory(result_partial)
#>   node time  S  I R
#> 1    2    3 99 NA 0
#> 2    4    3 99 NA 1
#> 3    2    5 98 NA 1
#> 4    4    5 99 NA 3

## Shortcut: If only 'time' and 'node' are specified, all
## compartments are recorded.
df <- data.frame(
  time = c(3, 5, 3, 5),
  node = c(2, 2, 4, 4)
)
punchcard(model) <- df

## Reset to record all data (equivalent to no template)
punchcard(model) <- NULL
result_reset <- run(model)
head(trajectory(result_reset))
#>   node time   S I R
#> 1    1    1 100 1 0
#> 2    2    1 100 3 0
#> 3    3    1 102 2 1
#> 4    4    1 103 4 0
#> 5    5    1 103 6 0
#> 6    6    1 105 6 0
```
