# Determine the number of compartments in a model

Extract the number of compartments from a `SimInf_model` object.
Compartments represent the distinct states an individual can occupy
within a node (e.g., Susceptible, Infected, Recovered). This count is
equivalent to the number of columns in the initial state vector `u0`.

## Usage

``` r
n_compartments(model)

# S4 method for class 'SimInf_model'
n_compartments(model)
```

## Arguments

- model:

  A `SimInf_model` object.

## Value

An integer scalar representing the total number of compartments in the
model.

## Examples

``` r
## Create an 'SIR' model with 3 compartments (S, I, R).
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

## Get the number of compartments.
n_compartments(model)
#> [1] 3
```
