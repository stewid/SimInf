# Set the shift matrix for a `SimInf_model` object

Utility function to set `events@N` in a `SimInf_model` object, see
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)

## Usage

``` r
shift_matrix(model) <- value

# S4 method for class 'SimInf_model'
shift_matrix(model) <- value
```

## Arguments

- model:

  The `model` to set the shift matrix `events@N`.

- value:

  A matrix.

## Value

`SimInf_model` object

## Examples

``` r
## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Set the shift matrix
shift_matrix(model) <- matrix(c(2, 1, 0), nrow = 3)

## Extract the shift matrix from the model
shift_matrix(model)
#>   1
#> S 2
#> I 1
#> R 0
```
