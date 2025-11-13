# Extract the shift matrix from a `SimInf_model` object

Utility function to extract the shift matrix `events@N` from a
`SimInf_model` object, see
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)

## Usage

``` r
shift_matrix(model)

# S4 method for class 'SimInf_model'
shift_matrix(model)
```

## Arguments

- model:

  The `model` to extract the shift matrix `events@N` from.

## Value

A mtrix.

## Examples

``` r
## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Extract the shift matrix from the model
shift_matrix(model)
#> <0 x 0 matrix>
```
