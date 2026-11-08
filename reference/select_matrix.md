# Extract the select matrix from a `SimInf_model` object

Get the select matrix `events@E` that determines which compartments an
event samples individuals from when processing *exit*, *enter*,
*internal transfer* and *external transfer* events. Each row in the
select matrix corresponds to one compartment in the model, and each
column defines a sampling pattern. The `select` value of an event
determines which column to use, and the values in the selected column
\\E\[, s\]\\ are used as weights when sampling individuals from the
compartments with non-zero entries. See
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for details on the event processing.

## Usage

``` r
select_matrix(model)

# S4 method for class 'SimInf_model'
select_matrix(model)
```

## Arguments

- model:

  The
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  object to extract the select matrix from.

## Value

The sparse select matrix `events@E` (of class
[`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html))
where each row corresponds to one compartment in the model.

## See also

[`select_matrix<-`](http://stewid.github.io/SimInf/reference/select_matrix-set.md)
to set the select matrix, and
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
