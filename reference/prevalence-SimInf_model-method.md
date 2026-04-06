# Calculate prevalence from a model object with trajectory data

Calculate the proportion of *cases* (specified on the left-hand side of
the formula) relative to the *population at risk* (specified on the
right-hand side) at different aggregation levels. The function supports
three levels of calculation:

- **Level 1 (Population Prevalence):** The proportion of cases in the
  total population at risk across all nodes.

- **Level 2 (Node Prevalence):** The proportion of nodes that contain at
  least one case, considering only nodes where the population at risk is
  greater than zero.

- **Level 3 (Within-Node Prevalence):** The proportion of cases relative
  to the population at risk calculated separately for each node.

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
## Create an 'SIR' model with 6 nodes.
u0 <- data.frame(S = 100:105, I = c(0, 1, 0, 2, 0, 3), R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Run the model. For reproducibility, we first call the
## set.seed() function and specify the number of threads to use
## since there is random sampling involved when picking individuals
## from the compartments.
set.seed(1)
set_num_threads(1)
result <- run(model)

## 1. Population Prevalence (level = 1, default)
## Proportion of infected individuals in the total population.
prevalence(result, I ~ S + I + R)
#>    time  prevalence
#> 1     1 0.009661836
#> 2     2 0.009661836
#> 3     3 0.008051530
#> 4     4 0.006441224
#> 5     5 0.004830918
#> 6     6 0.004830918
#> 7     7 0.003220612
#> 8     8 0.003220612
#> 9     9 0.004830918
#> 10   10 0.004830918

## Shorthand: '.' represents all compartments in the model (S + I + R).
prevalence(result, I ~ .)
#>    time  prevalence
#> 1     1 0.009661836
#> 2     2 0.009661836
#> 3     3 0.008051530
#> 4     4 0.006441224
#> 5     5 0.004830918
#> 6     6 0.004830918
#> 7     7 0.003220612
#> 8     8 0.003220612
#> 9     9 0.004830918
#> 10   10 0.004830918

## 2. Node Prevalence (level = 2)
## Proportion of nodes with at least one infected individual.
prevalence(result, I ~ S + I + R, level = 2)
#>    time prevalence
#> 1     1  0.5000000
#> 2     2  0.5000000
#> 3     3  0.5000000
#> 4     4  0.5000000
#> 5     5  0.5000000
#> 6     6  0.5000000
#> 7     7  0.3333333
#> 8     8  0.3333333
#> 9     9  0.3333333
#> 10   10  0.3333333

## 3. Within-Node Prevalence (level = 3)
## Prevalence calculated separately for each node.
prevalence(result, I ~ S + I + R, level = 3)
#>    node time  prevalence
#> 1     1    1 0.000000000
#> 2     2    1 0.009803922
#> 3     3    1 0.000000000
#> 4     4    1 0.019047619
#> 5     5    1 0.000000000
#> 6     6    1 0.027777778
#> 7     1    2 0.000000000
#> 8     2    2 0.009803922
#> 9     3    2 0.000000000
#> 10    4    2 0.019047619
#> 11    5    2 0.000000000
#> 12    6    2 0.027777778
#> 13    1    3 0.000000000
#> 14    2    3 0.009803922
#> 15    3    3 0.000000000
#> 16    4    3 0.019047619
#> 17    5    3 0.000000000
#> 18    6    3 0.018518519
#> 19    1    4 0.000000000
#> 20    2    4 0.009803922
#> 21    3    4 0.000000000
#> 22    4    4 0.019047619
#> 23    5    4 0.000000000
#> 24    6    4 0.009259259
#> 25    1    5 0.000000000
#> 26    2    5 0.009803922
#> 27    3    5 0.000000000
#> 28    4    5 0.009523810
#> 29    5    5 0.000000000
#> 30    6    5 0.009259259
#> 31    1    6 0.000000000
#> 32    2    6 0.009803922
#> 33    3    6 0.000000000
#> 34    4    6 0.009523810
#> 35    5    6 0.000000000
#> 36    6    6 0.009259259
#> 37    1    7 0.000000000
#> 38    2    7 0.009803922
#> 39    3    7 0.000000000
#> 40    4    7 0.009523810
#> 41    5    7 0.000000000
#> 42    6    7 0.000000000
#> 43    1    8 0.000000000
#> 44    2    8 0.009803922
#> 45    3    8 0.000000000
#> 46    4    8 0.009523810
#> 47    5    8 0.000000000
#> 48    6    8 0.000000000
#> 49    1    9 0.000000000
#> 50    2    9 0.009803922
#> 51    3    9 0.000000000
#> 52    4    9 0.019047619
#> 53    5    9 0.000000000
#> 54    6    9 0.000000000
#> 55    1   10 0.000000000
#> 56    2   10 0.009803922
#> 57    3   10 0.000000000
#> 58    4   10 0.019047619
#> 59    5   10 0.000000000
#> 60    6   10 0.000000000

## 4. Conditional Prevalence
## Calculate prevalence only in nodes where the number of
## recovered is zero.
prevalence(result, I ~ S + I + R | R == 0, level = 3)
#>    node time  prevalence
#> 1     1    1 0.000000000
#> 2     2    1 0.009803922
#> 3     3    1 0.000000000
#> 4     4    1 0.019047619
#> 5     5    1 0.000000000
#> 6     6    1 0.027777778
#> 7     1    2 0.000000000
#> 8     2    2 0.009803922
#> 9     3    2 0.000000000
#> 10    4    2 0.019047619
#> 11    5    2 0.000000000
#> 12    6    2         NaN
#> 13    1    3 0.000000000
#> 14    2    3 0.009803922
#> 15    3    3 0.000000000
#> 16    4    3 0.019047619
#> 17    5    3 0.000000000
#> 18    6    3         NaN
#> 19    1    4 0.000000000
#> 20    2    4 0.009803922
#> 21    3    4 0.000000000
#> 22    4    4         NaN
#> 23    5    4 0.000000000
#> 24    6    4         NaN
#> 25    1    5 0.000000000
#> 26    2    5 0.009803922
#> 27    3    5 0.000000000
#> 28    4    5         NaN
#> 29    5    5 0.000000000
#> 30    6    5         NaN
#> 31    1    6 0.000000000
#> 32    2    6 0.009803922
#> 33    3    6 0.000000000
#> 34    4    6         NaN
#> 35    5    6 0.000000000
#> 36    6    6         NaN
#> 37    1    7 0.000000000
#> 38    2    7 0.009803922
#> 39    3    7 0.000000000
#> 40    4    7         NaN
#> 41    5    7 0.000000000
#> 42    6    7         NaN
#> 43    1    8 0.000000000
#> 44    2    8 0.009803922
#> 45    3    8 0.000000000
#> 46    4    8         NaN
#> 47    5    8 0.000000000
#> 48    6    8         NaN
#> 49    1    9 0.000000000
#> 50    2    9 0.009803922
#> 51    3    9 0.000000000
#> 52    4    9         NaN
#> 53    5    9 0.000000000
#> 54    6    9         NaN
#> 55    1   10 0.000000000
#> 56    2   10 0.009803922
#> 57    3   10 0.000000000
#> 58    4   10         NaN
#> 59    5   10 0.000000000
#> 60    6   10         NaN
```
