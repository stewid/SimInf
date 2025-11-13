# Extract the select matrix from a `SimInf_model` object

Utility function to extract `events@E` from a `SimInf_model` object, see
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)

## Usage

``` r
select_matrix(model)

# S4 method for class 'SimInf_model'
select_matrix(model)
```

## Arguments

- model:

  The `model` to extract the select matrix `E` from.

## Value

[`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
object.

## Examples

``` r
## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Extract the select matrix from the model
select_matrix(model)
#> 3 x 4 sparse Matrix of class "dgCMatrix"
#>   1 2 3 4
#> S 1 . . 1
#> I . 1 . 1
#> R . . 1 1
```
