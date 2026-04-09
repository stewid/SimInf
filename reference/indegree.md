# Determine in-degree for each node in a model

Calculate the in-degree of each node based on **external transfer**
events (`"extTrans"`) in the model's schedules events.

## Usage

``` r
indegree(model)
```

## Arguments

- model:

  A `SimInf_model` object containing the event schedule.

## Value

An integer vector where each element corresponds to a node, containing
the count of unique source nodes sending individuals to it.

## Details

The in-degree is defined as the number of **unique source nodes** that
have sent individuals to the target node at least once. This metric
measures the connectivity of the network, indicating how many different
neighbors directly supply individuals to a specific node.

## Examples

``` r
## Create an 'SIR' model with example data.
model <- SIR(
  u0 = u0_SIR(),
  tspan = 1:1460,
  events = events_SIR(),
  beta = 0.16,
  gamma = 0.077
)

## Calculate the in-degree for each node.
deg <- indegree(model)

## View the in-degree for the first 6 nodes.
head(deg)
#> [1] 64 63 62 61 64 57

## Plot the distribution of in-degrees across all nodes.  This
## shows how many source nodes typically send individuals to a given
## node.
hist(
  deg,
  main = "Distribution of In-Degree (Unique Sources)",
  xlab = "Number of Unique Source Nodes"
)
```
