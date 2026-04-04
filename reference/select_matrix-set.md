# Set the select matrix for a `SimInf_model` object

Utility function to set `events@E` in a `SimInf_model` object, see
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md).

## Usage

``` r
select_matrix(model) <- value

# S4 method for class 'SimInf_model'
select_matrix(model) <- value
```

## Arguments

- model:

  The `SimInf_model` object to set the select matrix for.

- value:

  The new value for `E` in the model. `E` is a matrix to handle
  scheduled events, see
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md).
  Each row in `E` corresponds to one compartment in the model. The
  non-zero entries in a column indicate the compartments to include in
  an event. For the *exit*, *internal transfer* and *external transfer*
  events, the values in `E[, select]` are used as weights when sampling
  individuals without replacement, with probability proportional to the
  weight. For the *enter* event, the values in `E[, select]` are used as
  weights when determining which compartment to add individuals to. If
  the column `E[, select]` contains several non-zero entries, the
  compartment is sampled with probability proportional to the weight in
  `E[, select]`. The select matrix `E` can either be specified as a
  `matrix`, or as a `data.frame`. When `E` is specified as a
  `data.frame`, it must have one column named `compartment` that defines
  which compartment is referred to, and one column `select` that defines
  the column in `E`. In addition, the `data.frame` can contain an
  optional column named `value` with the value in `E`. When the `value`
  column is missing, `1` is used as the default value.

## Examples

``` r
## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Set the select matrix
select_matrix(model) <- matrix(c(1, 0, 0, 1, 1, 1, 0, 0, 1), nrow = 3)

## Extract the select matrix from the model
select_matrix(model)
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>   1 2 3
#> S 1 1 .
#> I . 1 .
#> R . 1 1
```
