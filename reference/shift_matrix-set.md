# Set the shift matrix for a `SimInf_model` object

Utility function to set `events@N` in a `SimInf_model` object, see
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md).

## Usage

``` r
shift_matrix(model) <- value

# S4 method for class 'SimInf_model'
shift_matrix(model) <- value
```

## Arguments

- model:

  The `SimInf_model` object to set the shift matrix for.

- value:

  The new value for `N` in the model. `N` is a matrix to handle
  scheduled events, see
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md).
  Each row in `N` corresponds to one compartment in the model. The
  values in a column define how to move sampled individuals before
  adding them to the destination. Let `q <- shift`, then each non-zero
  entry in `N[, q]` defines the number of rows to move sampled
  individuals from that compartment i.e., sampled individuals from
  compartment `p` are moved to compartment `N[p, q] + p`, where
  `1 <= N[p, q] + p <= N_compartments`. This matrix is used for *enter*,
  *internal transfer* and *external transfer* events. The shift matrix
  `N` can either be specified as a `matrix`, or as a `data.frame`. When
  `N` is specified as a `data.frame`, it must have one column named
  `compartment` that defines which compartment is referred to, and one
  column `shift` that defines the column in `N`. In addition, the
  `data.frame` must contain a column named `value` with the integer
  value in `N`.

## Value

`SimInf_model` object

## Examples

``` r
## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Set the shift matrix.
shift_matrix(model) <- matrix(c(2, 1, 0), nrow = 3)

## Extract the shift matrix from the model.
shift_matrix(model)
#>   1
#> S 2
#> I 1
#> R 0

## Set the shift matrix using a data.frame instead.
shift_matrix(model) <- data.frame(
    compartment = c("S", "I", "R"),
    shift = c(1, 1, 1),
    value = c(2, 1, 0))

## Extract the shift matrix from the model.
shift_matrix(model)
#>   1
#> S 2
#> I 1
#> R 0
```
