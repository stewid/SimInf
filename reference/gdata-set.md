# Set a global data parameter for a `SimInf_model` object

The global data is a numeric vector that is common to all nodes. The
global data vector is passed as an argument to the transition rate
functions and the post time step function.

## Usage

``` r
gdata(model, parameter) <- value

# S4 method for class 'SimInf_model'
gdata(model, parameter) <- value
```

## Arguments

- model:

  The `model` to set a global model parameter for.

- parameter:

  The name of the parameter to set.

- value:

  A numeric value.

## Value

a `SimInf_model` object

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
