# Determine the number of nodes in a model

Determine the number of nodes in a model

## Usage

``` r
n_nodes(model)

# S4 method for class 'SimInf_model'
n_nodes(model)

# S4 method for class 'SimInf_pfilter'
n_nodes(model)

# S4 method for class 'SimInf_pmcmc'
n_nodes(model)
```

## Arguments

- model:

  the `model` object to extract the number of nodes from.

## Value

the number of nodes in the model.

## Examples

``` r
## Create an 'SIR' model with 100 nodes, with 99 susceptible,
## 1 infected and 0 recovered in each node.
u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Display the number of nodes in the model.
n_nodes(model)
#> [1] 100
```
