# Determine out-degree for each node in a model

Calculate the out-degree of each node based on **external transfer**
events (`"extTrans"`) in the model's scheduled events.

## Usage

``` r
outdegree(model)
```

## Arguments

- model:

  A `SimInf_model` object containing the scheduled events.

## Value

An integer vector where each element corresponds to a node, containing
the count of unique destination nodes receiving individuals from it.

## Details

The out-degree is defined as the number of **unique destination nodes**
that receive individuals from the source node at least once. This metric
measures the connectivity of the network, indicating how many different
neighbors a specific node directly sends individuals to.

## See also

[`indegree`](http://stewid.github.io/SimInf/reference/indegree.md) for
calculating the number of unique source nodes that send individuals to a
node.
[`events_SIR`](http://stewid.github.io/SimInf/reference/events_SIR.md)
for example event data used in network analysis.
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for details on the structure of scheduled events.

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

## Calculate the out-degree for each node.
deg <- outdegree(model)

## View the out-degree for the first 6 nodes.
head(deg)
#> [1] 61 60 64 64 67 66

## Plot the distribution of out-degrees across all nodes.
## This shows how many destination nodes typically receive
## individuals from a given node.
hist(
  deg,
  main = "Distribution of Out-Degree (Unique Destinations)",
  xlab = "Number of Unique Destination Nodes"
)
```
