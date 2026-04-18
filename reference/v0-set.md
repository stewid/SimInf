# Update the initial continuous state (`v0`) in each node

Replace the initial continuous state vector (`v0`) of a `SimInf_model`
object with new data. This allows you to modify the starting conditions
of continuous variables (e.g., environmental pathogen concentration)
without recreating the object.

## Usage

``` r
v0(model) <- value

# S4 method for class 'SimInf_model'
v0(model) <- value
```

## Arguments

- model:

  A `SimInf_model` object.

- value:

  An object containing the new initial continuous state. Can be a
  `data.frame`, `matrix`, or `named numeric vector`. Non-data.frame
  inputs will be coerced to a `data.frame`.

## Value

The modified `SimInf_model` object.

## Details

The `value` argument accepts a `data.frame`, `matrix`, or
`named numeric vector`. If the input is not a `data.frame`, it will be
automatically coerced to one. The function handles the following
formats:

- **Single Node**: If `value` is a named vector or a one-row
  matrix/data.frame, it is applied to the single node in the model.

- **Multiple Nodes**: If `value` is a matrix or data.frame with multiple
  rows, each row corresponds to one node. The number of rows must
  exactly match the number of nodes in the `model`.

- **Column Matching**: Column names must match the continuous state
  variable names defined in the model (e.g., `"phi"`). Only matching
  columns are used; extra columns are ignored, and missing variables
  will trigger an error.

The function validates the input and ensures the new state is consistent
with the model structure before updating.

## See also

`u0<-` for updating the initial discrete compartment state.

## Examples

``` r
## For reproducibility, set the seed.
set.seed(22)

## Create an 'SISe' model with no infected individuals and no
## infectious pressure (phi = 0).
model <- SISe(
  u0 = data.frame(S = 100, I = 0),
  tspan = 1:100,
  phi = 0,
  upsilon = 0.02,
  gamma = 0.1,
  alpha = 1,
  epsilon = 0,
  beta_t1 = 0.15,
  beta_t2 = 0.15,
  beta_t3 = 0.15,
  beta_t4 = 0.15,
  end_t1 = 91,
  end_t2 = 182,
  end_t3 = 273,
  end_t4 = 365
)

## Run the 'SISe' model and plot the result.
result <- run(model)
plot(result)


## Update the infectious pressure 'phi' in 'v0' using a named
## vector.  (Automatically coerced to one row for the single
## node).
v0(model) <- c(phi = 1)
result <- run(model)
plot(result)


## For a multi-node model, use a data.frame with multiple rows:
## v0(model) <- data.frame(phi = c(1.0, 0.5, 0.0))
```
