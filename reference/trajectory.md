# Generic function to extract data from a simulated trajectory

Generic function to extract data from a simulated trajectory

## Usage

``` r
trajectory(model, compartments = NULL, index = NULL, ...)
```

## Arguments

- model:

  the object to extract the trajectory from.

- compartments:

  specify the names of the compartments to extract data from. The
  compartments can be specified as a character vector e.g.
  `compartments = c('S', 'I', 'R')`, or as a formula e.g.
  `compartments = ~S+I+R` (see ‘Examples’). Default
  (`compartments=NULL`) is to extract the number of individuals in each
  compartment i.e. the data from all discrete state compartments in the
  model. In models that also have continuous state variables e.g. the
  `SISe` model, they are also included.

- index:

  indices specifying the subset of nodes to include when extracting
  data. Default (`index = NULL`) is to extract data from all nodes.

- ...:

  Additional arguments, see
  [`trajectory,SimInf_model-method`](http://stewid.github.io/SimInf/reference/trajectory-SimInf_model-method.md)
