# Determine out-degree for each node in a model

The number nodes that are connected with *external transfer* events from
each node.

## Usage

``` r
outdegree(model)
```

## Arguments

- model:

  determine out-degree for each node in the model.

## Value

vector with out-degree for each node.

## Examples

``` r
## Create an 'SIR' model with 1600 nodes and initialize
## it with example data.
model <- SIR(u0 = u0_SIR(), tspan = 1:1460, events = events_SIR(),
             beta   = 0.16, gamma  = 0.077)

## Display outdegree for each node in the model.
plot(outdegree(model))
```
