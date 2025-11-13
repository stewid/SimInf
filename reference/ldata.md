# Extract local data from a node

The local data is a numeric vector that is specific to a node. The local
data vector is passed as an argument to the transition rate functions
and the post time step function.

## Usage

``` r
ldata(model, node)

# S4 method for class 'SimInf_model'
ldata(model, node)
```

## Arguments

- model:

  The `model` to get local data from.

- node:

  index to node to extract local data from.

## Value

a numeric vector

## Examples

``` r
## Create an 'SISe' model with 1600 nodes.
model <- SISe(u0 = u0_SISe(), tspan = 1:100, events = events_SISe(),
              phi = 0, upsilon = 1.8e-2, gamma = 0.1, alpha = 1,
              beta_t1 = 1.0e-1, beta_t2 = 1.0e-1, beta_t3 = 1.25e-1,
              beta_t4 = 1.25e-1, end_t1 = c(91, 101), end_t2 = c(182, 185),
              end_t3 = c(273, 275), end_t4 = c(365, 360), epsilon = 0)

## Display local data from the first two nodes.
ldata(model, node = 1)
#> end_t1 end_t2 end_t3 end_t4 
#>     91    182    273    365 
ldata(model, node = 2)
#> end_t1 end_t2 end_t3 end_t4 
#>    101    185    275    360 
```
