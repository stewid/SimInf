# Determine the number of nodes in a model

Extract the number of nodes from a `SimInf_model` object. A node
represents a distinct sub-population in the model, but its definition is
determined by the modeller. For example, a node can represent a cattle
herd, a pen within a herd, or even a single individual, depending on the
research question and the scale of the study. This count is equivalent
to the number of rows in the initial state vector `u0`.

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

  A `SimInf_model` object.

## Value

An integer scalar representing the total number of nodes in the model.

## Examples

``` r
## Create an 'SIR' model with 100 nodes.
u0 <- data.frame(
  S = rep(99, 100),
  I = rep(1, 100),
  R = rep(0, 100)
)

model <- SIR(
  u0 = u0,
  tspan = 1:10,
  beta = 0.16,
  gamma = 0.077
)

## Get the number of nodes.
n_nodes(model)
#> [1] 100
```
