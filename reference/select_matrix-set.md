# Set the select matrix for a `SimInf_model` object

Utility function to set `events@E` in a `SimInf_model` object, see
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)

## Usage

``` r
select_matrix(model) <- value

# S4 method for class 'SimInf_model'
select_matrix(model) <- value
```

## Arguments

- model:

  The `model` to set the select matrix for.

- value:

  A matrix.

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
