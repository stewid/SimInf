# Extract filtered trajectories from fitting a PMCMC algorithm

Extract filtered trajectories from a particle Markov chain Monte Carlo
algorithm.

## Usage

``` r
# S4 method for class 'SimInf_pmcmc'
trajectory(model, compartments, index, start = 1, end = NULL, thin = 1)
```

## Arguments

- model:

  the `SimInf_pmcmc` object to extract the filtered trajectories from.

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

- start:

  The start iteration to remove some burn-in iterations. Default is
  `start = 1`.

- end:

  the last iteration to include. Default is `NULL` which set `end` to
  the last iteration in the chain.

- thin:

  keep every `thin` iteration after the `start` iteration. Default is
  `thin = 1`, i.e., keep every iteration.

## Value

A `data.frame` where the first column is the `iteration` and the
remaining columns are the result from calling
[`trajectory,SimInf_model-method`](http://stewid.github.io/SimInf/reference/trajectory-SimInf_model-method.md)
with the arguments `compartments` and `index` for each iteration.
