# Determine in-degree for each node in a model

The number of nodes with inward *external transfer* events to each node.

## Usage

``` r
indegree(model)
```

## Arguments

- model:

  determine in-degree for each node in the model.

## Value

vector with in-degree for each node.

## Examples

``` r
## Create an 'SIR' model with 1600 nodes and initialize
## it with example data.
model <- SIR(u0 = u0_SIR(), tspan = 1:1460, events = events_SIR(),
             beta   = 0.16, gamma  = 0.077)

## Display indegree for each node in the model.
plot(indegree(model))
```
