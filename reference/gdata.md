# Extract global data from a `SimInf_model` object

The global data is a numeric vector that is common to all nodes. The
global data vector is passed as an argument to the transition rate
functions and the post time step function.

## Usage

``` r
gdata(model)

# S4 method for class 'SimInf_model'
gdata(model)
```

## Arguments

- model:

  The `model` to get global data from.

## Value

a numeric vector

## Examples

``` r
## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Set 'beta' to a new value
gdata(model, "beta") <- 2

## Extract the global data vector that is common to all nodes
gdata(model)
#> beta 
#>    2 
```
