# Set the select matrix for a `SimInf_model` object

Replace the select matrix `events@E` that determines which compartments
an event samples individuals from when processing *exit*, *enter*,
*internal transfer* and *external transfer* events. The values in the
selected column are used as weights when sampling individuals from the
compartments with non-zero entries. The select matrix can be specified
as a `matrix`, or as a `data.frame` with the columns `compartment` and
`select`, and an optional column `value` (default `1`); see ‘Arguments’
for details, and
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for how the select matrix is used when processing events.

## Usage

``` r
select_matrix(model) <- value

# S4 method for class 'SimInf_model'
select_matrix(model) <- value
```

## Arguments

- model:

  The
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  object to set the select matrix for.

- value:

  The new value for `E` in the model. `E` is a matrix to handle
  scheduled events, see
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md).
  Each row in `E` corresponds to one compartment in the model. The
  values in a column are used as weights when sampling individuals from
  the compartments with non-zero entries when processing an event. Let
  `s <- select`, then the non-zero entries in `E[, s]` determine which
  compartments to sample individuals from. This matrix is used for
  *exit*, *enter*, *internal transfer* and *external transfer* events.
  The select matrix `E` can either be specified as a `matrix`, or as a
  `data.frame`. When `E` is specified as a `data.frame`, it must have
  one column named `compartment` that defines which compartment is
  referred to, and one column `select` that defines the column in `E`.
  In addition, the `data.frame` may contain a column named `value` with
  the weight in `E` (`1` if the column is omitted). Note that the
  `data.frame` replaces the entire select matrix: entries not included
  are set to zero, so all columns that the scheduled events refer to
  must be specified.

## Value

The
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object with the updated select matrix.

## See also

[`select_matrix`](http://stewid.github.io/SimInf/reference/select_matrix.md)
to extract the select matrix, and
[`shift_matrix`](http://stewid.github.io/SimInf/reference/shift_matrix.md)
to extract the shift matrix. The vignette "Scheduled events" for a
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

## Extract the select matrix from the model.
select_matrix(model)
#> 3 x 4 sparse Matrix of class "dgCMatrix"
#>   1 2 3 4
#> S 1 . . 1
#> I . 1 . 1
#> R . . 1 1

## Change the select matrix so that individuals in the 'I'
## compartment are twice as likely to be selected when
## processing an event with 'select' = 4. Note that the
## data.frame replaces the entire select matrix, so the entries
## in columns 1, 2 and 3 must also be specified.
select_matrix(model) <- data.frame(
  compartment = c("S", "I", "R", "S", "I", "R"),
  select = c(1, 2, 3, 4, 4, 4),
  value = c(1, 1, 1, 1, 2, 1)
)

## Extract the updated select matrix from the model.
select_matrix(model)
#> 3 x 4 sparse Matrix of class "dgCMatrix"
#>   1 2 3 4
#> S 1 . . 1
#> I . 1 . 2
#> R . . 1 1
```
