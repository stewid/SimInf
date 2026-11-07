# Extract the shift matrix from a `SimInf_model` object

Get the shift matrix `events@N` that determines how individuals are
shifted between compartments when processing *enter*, *internal
transfer* and *external transfer* events. Each row in the shift matrix
corresponds to one compartment in the model, and each column defines a
transfer pattern. The `shift` value of an event selects which column to
use, and the value \\N\[p, q\]\\ indicates how many compartments to move
(or place, for enter events) individuals from compartment \\p\\. See
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for details on the event processing.

## Usage

``` r
shift_matrix(model)

# S4 method for class 'SimInf_model'
shift_matrix(model)
```

## Arguments

- model:

  The
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  object to extract the shift matrix from.

## Value

The integer shift matrix `events@N` where each row corresponds to one
compartment in the model.

## See also

[`shift_matrix<-`](http://stewid.github.io/SimInf/reference/shift_matrix-set.md)
to set the shift matrix, and
[`select_matrix`](http://stewid.github.io/SimInf/reference/select_matrix.md)
to extract the select matrix. The vignette "Scheduled events" for a
tutorial on scheduled events, see
`vignette("Scheduled events", package = "SimInf")`.

## Examples

``` r
## Create an 'SIR' model.
model <- SIR(
  u0 = data.frame(
    S = 100,
    I = 10,
    R = 0),
  tspan = 1:5,
  beta = 0.16,
  gamma = 0.077
)

## Extract the shift matrix from the model.
shift_matrix(model)
#> <0 x 0 matrix>

## Change the shift matrix so that susceptible individuals are
## moved to the 'R' compartment (compartment 1 + 2 = 3) when
## processing an event with 'shift' = 1.
shift_matrix(model) <- data.frame(
  compartment = "S",
  shift = 1,
  value = 2
)

## Extract the updated shift matrix from the model.
shift_matrix(model)
#>   1
#> S 2
#> I 0
#> R 0
```
