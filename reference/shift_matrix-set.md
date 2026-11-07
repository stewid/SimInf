# Set the shift matrix for a `SimInf_model` object

Replace the shift matrix `events@N` that determines how individuals are
shifted between compartments when processing *enter*, *internal
transfer* and *external transfer* events. The shift matrix can be
specified as a `matrix`, or as a `data.frame` with the columns
`compartment`, `shift`, and `value`; see ‘Arguments’ for details, and
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for how the shift matrix is used when processing events.

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
  entry in `N[, q]` defines the number of compartments to move sampled
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

## See also

[`shift_matrix`](http://stewid.github.io/SimInf/reference/shift_matrix.md)
to extract the shift matrix, and
[`select_matrix`](http://stewid.github.io/SimInf/reference/select_matrix.md)
to extract the select matrix. The vignette "Scheduled events" for a
tutorial on scheduled events, see
`vignette("Scheduled events", package = "SimInf")`.

## Examples

``` r
## Create an 'SIR' model and turn off the disease dynamics to
## focus on the scheduled event.
model <- SIR(
  u0 = data.frame(
    S = 100,
    I = 10,
    R = 0),
  tspan = 1:20,
  beta = 0,
  gamma = 0,
  events = data.frame(
    event = "intTrans",
    time = 10,
    node = 1,
    dest = 0,
    n = 30,
    proportion = 0,
    select = 1,
    shift = 1
  )
)

## Vaccinate 30 susceptible individuals at t = 10 by moving
## them to the 'R' compartment (compartment 1 + 2 = 3).
shift_matrix(model) <- data.frame(
  compartment = "S",
  shift = 1,
  value = 2
)

## Run the model and extract the trajectory.
result <- run(model)
trajectory(result)
#>    node time   S  I  R
#> 1     1    1 100 10  0
#> 2     1    2 100 10  0
#> 3     1    3 100 10  0
#> 4     1    4 100 10  0
#> 5     1    5 100 10  0
#> 6     1    6 100 10  0
#> 7     1    7 100 10  0
#> 8     1    8 100 10  0
#> 9     1    9 100 10  0
#> 10    1   10  70 10 30
#> 11    1   11  70 10 30
#> 12    1   12  70 10 30
#> 13    1   13  70 10 30
#> 14    1   14  70 10 30
#> 15    1   15  70 10 30
#> 16    1   16  70 10 30
#> 17    1   17  70 10 30
#> 18    1   18  70 10 30
#> 19    1   19  70 10 30
#> 20    1   20  70 10 30
```
