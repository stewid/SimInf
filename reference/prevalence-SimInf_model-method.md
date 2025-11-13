# Calculate prevalence from a model object with trajectory data

Calculate the proportion of individuals with disease in the population,
or the proportion of nodes with at least one diseased individual, or the
proportion of individuals with disease in each node.

## Usage

``` r
# S4 method for class 'SimInf_model'
prevalence(model, formula, level, index, format = c("data.frame", "matrix"))
```

## Arguments

- model:

  The `model` with trajectory data to calculate the prevalence from.

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

## Examples

``` r
## Create an 'SIR' model with 6 nodes and initialize
## it to run over 10 days.
u0 <- data.frame(S = 100:105, I = c(0, 1, 0, 2, 0, 3), R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Determine the proportion of infected individuals (cases)
## in the population at the time-points in 'tspan'.
prevalence(result, I ~ S + I + R)
#>    time prevalence
#> 1     1 0.00805153
#> 2     2 0.01127214
#> 3     3 0.01288245
#> 4     4 0.01449275
#> 5     5 0.01610306
#> 6     6 0.01771337
#> 7     7 0.01771337
#> 8     8 0.02415459
#> 9     9 0.02576490
#> 10   10 0.02737520

## Identical result is obtained with the shorthand 'I~.'
prevalence(result, I ~ .)
#>    time prevalence
#> 1     1 0.00805153
#> 2     2 0.01127214
#> 3     3 0.01288245
#> 4     4 0.01449275
#> 5     5 0.01610306
#> 6     6 0.01771337
#> 7     7 0.01771337
#> 8     8 0.02415459
#> 9     9 0.02576490
#> 10   10 0.02737520

## Determine the proportion of nodes with infected individuals at
## the time-points in 'tspan'.
prevalence(result, I ~ S + I + R, level = 2)
#>    time prevalence
#> 1     1  0.3333333
#> 2     2  0.3333333
#> 3     3  0.3333333
#> 4     4  0.3333333
#> 5     5  0.3333333
#> 6     6  0.3333333
#> 7     7  0.3333333
#> 8     8  0.3333333
#> 9     9  0.3333333
#> 10   10  0.3333333

## Determine the proportion of infected individuals in each node
## at the time-points in 'tspan'.
prevalence(result, I ~ S + I + R, level = 3)
#>    node time prevalence
#> 1     1    1 0.00000000
#> 2     2    1 0.00000000
#> 3     3    1 0.00000000
#> 4     4    1 0.01904762
#> 5     5    1 0.00000000
#> 6     6    1 0.02777778
#> 7     1    2 0.00000000
#> 8     2    2 0.00000000
#> 9     3    2 0.00000000
#> 10    4    2 0.01904762
#> 11    5    2 0.00000000
#> 12    6    2 0.04629630
#> 13    1    3 0.00000000
#> 14    2    3 0.00000000
#> 15    3    3 0.00000000
#> 16    4    3 0.01904762
#> 17    5    3 0.00000000
#> 18    6    3 0.05555556
#> 19    1    4 0.00000000
#> 20    2    4 0.00000000
#> 21    3    4 0.00000000
#> 22    4    4 0.02857143
#> 23    5    4 0.00000000
#> 24    6    4 0.05555556
#> 25    1    5 0.00000000
#> 26    2    5 0.00000000
#> 27    3    5 0.00000000
#> 28    4    5 0.02857143
#> 29    5    5 0.00000000
#> 30    6    5 0.06481481
#> 31    1    6 0.00000000
#> 32    2    6 0.00000000
#> 33    3    6 0.00000000
#> 34    4    6 0.03809524
#> 35    5    6 0.00000000
#> 36    6    6 0.06481481
#> 37    1    7 0.00000000
#> 38    2    7 0.00000000
#> 39    3    7 0.00000000
#> 40    4    7 0.04761905
#> 41    5    7 0.00000000
#> 42    6    7 0.05555556
#> 43    1    8 0.00000000
#> 44    2    8 0.00000000
#> 45    3    8 0.00000000
#> 46    4    8 0.05714286
#> 47    5    8 0.00000000
#> 48    6    8 0.08333333
#> 49    1    9 0.00000000
#> 50    2    9 0.00000000
#> 51    3    9 0.00000000
#> 52    4    9 0.05714286
#> 53    5    9 0.00000000
#> 54    6    9 0.09259259
#> 55    1   10 0.00000000
#> 56    2   10 0.00000000
#> 57    3   10 0.00000000
#> 58    4   10 0.04761905
#> 59    5   10 0.00000000
#> 60    6   10 0.11111111

## Determine the proportion of infected individuals in each node
## at the time-points in 'tspan' when the number of recovered is
## zero.
prevalence(result, I ~ S + I + R | R == 0, level = 3)
#>    node time prevalence
#> 1     1    1 0.00000000
#> 2     2    1        NaN
#> 3     3    1 0.00000000
#> 4     4    1 0.01904762
#> 5     5    1 0.00000000
#> 6     6    1 0.02777778
#> 7     1    2 0.00000000
#> 8     2    2        NaN
#> 9     3    2 0.00000000
#> 10    4    2 0.01904762
#> 11    5    2 0.00000000
#> 12    6    2 0.04629630
#> 13    1    3 0.00000000
#> 14    2    3        NaN
#> 15    3    3 0.00000000
#> 16    4    3 0.01904762
#> 17    5    3 0.00000000
#> 18    6    3 0.05555556
#> 19    1    4 0.00000000
#> 20    2    4        NaN
#> 21    3    4 0.00000000
#> 22    4    4 0.02857143
#> 23    5    4 0.00000000
#> 24    6    4        NaN
#> 25    1    5 0.00000000
#> 26    2    5        NaN
#> 27    3    5 0.00000000
#> 28    4    5 0.02857143
#> 29    5    5 0.00000000
#> 30    6    5        NaN
#> 31    1    6 0.00000000
#> 32    2    6        NaN
#> 33    3    6 0.00000000
#> 34    4    6 0.03809524
#> 35    5    6 0.00000000
#> 36    6    6        NaN
#> 37    1    7 0.00000000
#> 38    2    7        NaN
#> 39    3    7 0.00000000
#> 40    4    7 0.04761905
#> 41    5    7 0.00000000
#> 42    6    7        NaN
#> 43    1    8 0.00000000
#> 44    2    8        NaN
#> 45    3    8 0.00000000
#> 46    4    8 0.05714286
#> 47    5    8 0.00000000
#> 48    6    8        NaN
#> 49    1    9 0.00000000
#> 50    2    9        NaN
#> 51    3    9 0.00000000
#> 52    4    9        NaN
#> 53    5    9 0.00000000
#> 54    6    9        NaN
#> 55    1   10 0.00000000
#> 56    2   10        NaN
#> 57    3   10 0.00000000
#> 58    4   10        NaN
#> 59    5   10 0.00000000
#> 60    6   10        NaN
```
