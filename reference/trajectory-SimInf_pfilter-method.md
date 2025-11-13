# Extract filtered trajectory from running a particle filter

Extract filtered trajectory from running a particle filter

## Usage

``` r
# S4 method for class 'SimInf_pfilter'
trajectory(model, compartments, index, format = c("data.frame", "matrix"))
```

## Arguments

- model:

  the `SimInf_pfilter` object to extract the result from.

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

- format:

  the default (`format = "data.frame"`) is to generate a `data.frame`
  with one row per node and time-step with the number of individuals in
  each compartment. Using `format = "matrix"` returns the result as a
  matrix, which is the internal format (see ‘Details’ in
  [`trajectory,SimInf_model-method`](http://stewid.github.io/SimInf/reference/trajectory-SimInf_model-method.md)).

## Value

A `data.frame` if `format = "data.frame"`, else a matrix.
