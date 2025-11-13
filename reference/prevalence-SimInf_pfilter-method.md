# Extract prevalence from running a particle filter

Extract prevalence from running a particle filter

## Usage

``` r
# S4 method for class 'SimInf_pfilter'
prevalence(model, formula, level, index, format = c("data.frame", "matrix"))
```

## Arguments

- model:

  the `SimInf_pfilter` object to extract the prevalence from.

- formula:

  A formula that specifies the compartments that define the cases with a
  disease or that have a specific characteristic (numerator), and the
  compartments that define the entire population of interest
  (denominator). The left-hand-side of the formula defines the cases,
  and the right-hand-side defines the population, for example, `I~S+I+R`
  in a ‘SIR’ model (see ‘Examples’). The `.` (dot) is expanded to all
  compartments, for example, `I~.` is expanded to `I~S+I+R` in a ‘SIR’
  model (see ‘Examples’). The formula can also contain a condition
  (indicated by `|`) for each node and time step to further control the
  population to include in the calculation, for example,
  `I ~ . | R == 0` to calculate the prevalence when the recovered is
  zero in a ‘SIR’ model. The condition must evaluate to `TRUE` or
  `FALSE` in each node and time step. Please note, if the denominator is
  zero, the prevalence is `NaN`. Additionally, when `level=3`
  (within-node prevalence) and the formula contains a condition that
  evaluates to `FALSE`, the prevalence is also `NaN`.

- level:

  The level at which the prevalence is calculated at each time point in
  `tspan`. 1 (population prevalence): calculates the proportion of the
  individuals (cases) in the population. 2 (node prevalence): calculates
  the proportion of nodes with at least one case. 3 (within-node
  prevalence): calculates the proportion of cases within each node.
  Default is `1`.

- index:

  indices specifying the subset of nodes to include when extracting
  data. Default (`index = NULL`) is to extract data from all nodes.

- format:

  The default (`format = "data.frame"`) is to generate a `data.frame`
  with one row per time-step with the prevalence. Using
  `format = "matrix"` returns the result as a matrix.

## Value

A `data.frame` if `format = "data.frame"`, else a matrix.
