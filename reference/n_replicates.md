# Determine the number of replicates in a model

Extract the number of replicates from a `SimInf_model` object.
Replicates are independent copies of the model state used by the
**particle filter**
([`pfilter`](http://stewid.github.io/SimInf/reference/pfilter.md)). This
value is set **internally** by `pfilter` and is not specified when
creating a standard model. If the model has not been processed by a
filter, the value will be 1.

## Usage

``` r
n_replicates(model)

# S4 method for class 'SimInf_model'
n_replicates(model)
```

## Arguments

- model:

  A `SimInf_model` object.

## Value

An integer scalar representing the number of replicates in the model.

## Examples

``` r
## Create a standard 'SIR' model.
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

## Get the number of replicates (default is 1 for standard
## models).
n_replicates(model)
#> [1] 1
```
