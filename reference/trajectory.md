# Extract data from a simulated trajectory

Extract the number of individuals in each compartment in every node
after generating a single stochastic trajectory with
[`run`](http://stewid.github.io/SimInf/reference/run.md). For a
[`SimInf_pfilter`](http://stewid.github.io/SimInf/reference/SimInf_pfilter-class.md)
object, extract the filtered trajectory from running a particle filter.

## Usage

``` r
trajectory(
  model,
  compartments = NULL,
  index = NULL,
  format = c("data.frame", "matrix")
)

# S4 method for class 'SimInf_model'
trajectory(
  model,
  compartments = NULL,
  index = NULL,
  format = c("data.frame", "matrix")
)

# S4 method for class 'SimInf_pfilter'
trajectory(
  model,
  compartments = NULL,
  index = NULL,
  format = c("data.frame", "matrix")
)
```

## Arguments

- model:

  the `SimInf_model` object to extract the result from.

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
  each compartment. When the model contains multiple replicates of each
  node, the `data.frame` also contains one column `replicate`. Using
  `format = "matrix"` returns the result as a matrix, which is the
  internal format (see ‘Details’).

## Value

A `data.frame` if `format = "data.frame"`, else a matrix.

## Internal format of the discrete state variables

When `format = "matrix"`, the result is the internal matrix `U`. The
column `U[, j]` contains the number of individuals in each compartment
at `tspan[j]`, where the first \\N_c\\ rows contain node 1, the next
\\N_c\\ rows contain node 2, and so on, where \\N_c\\ is the number of
compartments in the model. Consequently, the dimension of the matrix is
\\N_n N_c \times\\ `length(tspan)`, where \\N_n\\ is the number of
nodes. The internal format also allows replicates of each node. Each
replicate adds new columns to `U`, so that the data for each replicate
is in blocks of `length(tspan)` columns.

## Internal format of the continuous state variables

When `format = "matrix"`, the result is the internal matrix `V`. The
column `V[, j]` contains the real-valued continuous state of the system
at `tspan[j]`, where the rows are grouped per node and the number of
rows per node is determined by the number of continuous state variables
in the model. Consequently, the dimension of the matrix is \\N_n N_d
\times\\ `length(tspan)`, where \\N_n\\ is the number of nodes and
\\N_d\\ is the number of continuous state variables. The internal format
also allows replicates of each node. Each replicate adds new columns to
`V`, so that the data for each replicate is in blocks of `length(tspan)`
columns.

## See also

[`run`](http://stewid.github.io/SimInf/reference/run.md) for generating
a trajectory, and
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence-SimInf_model-method.md)
for calculating prevalence from trajectory data. The vignette
"Post-process data in a trajectory" for a tutorial on post-processing,
see `vignette("Post-process data in a trajectory", package = "SimInf")`.

## Examples

``` r
## For reproducibility, specify the number of threads.
set_num_threads(1)

## Create an 'SIR' model with 6 nodes and initialize
## it to run over 10 days.
model <- SIR(
  u0 = data.frame(
    S = 100:105,
    I = 1:6,
    R = rep(0, 6)),
  tspan = 1:10,
  beta = 0.16,
  gamma = 0.077
)

## Run the model with a fixed seed for reproducibility.
result <- run(model, seed = 123)

## Extract the number of individuals in each compartment at the
## time-points in 'tspan'.
trajectory(result)
#>    node time   S  I  R
#> 1     1    1 100  1  0
#> 2     2    1 101  2  0
#> 3     3    1 102  3  0
#> 4     4    1 103  4  0
#> 5     5    1 100  9  0
#> 6     6    1 104  7  0
#> 7     1    2 100  1  0
#> 8     2    2 101  2  0
#> 9     3    2 102  2  1
#> 10    4    2 103  3  1
#> 11    5    2  99  9  1
#> 12    6    2 104  6  1
#> 13    1    3 100  0  1
#> 14    2    3 101  2  0
#> 15    3    3 102  2  1
#> 16    4    3 103  3  1
#> 17    5    3  97  8  4
#> 18    6    3 104  6  1
#> 19    1    4 100  0  1
#> 20    2    4 101  2  0
#> 21    3    4 102  1  2
#> 22    4    4 103  3  1
#> 23    5    4  95 10  4
#> 24    6    4 103  7  1
#> 25    1    5 100  0  1
#> 26    2    5 101  2  0
#> 27    3    5 102  1  2
#> 28    4    5 101  5  1
#> 29    5    5  94 11  4
#> 30    6    5 103  7  1
#> 31    1    6 100  0  1
#> 32    2    6 101  2  0
#> 33    3    6 102  0  3
#> 34    4    6 100  6  1
#> 35    5    6  92 11  6
#> 36    6    6 102  7  2
#> 37    1    7 100  0  1
#> 38    2    7 101  1  1
#> 39    3    7 102  0  3
#> 40    4    7 100  5  2
#> 41    5    7  89 13  7
#> 42    6    7 102  5  4
#> 43    1    8 100  0  1
#> 44    2    8 101  1  1
#> 45    3    8 102  0  3
#> 46    4    8  98  7  2
#> 47    5    8  88 14  7
#> 48    6    8 102  5  4
#> 49    1    9 100  0  1
#> 50    2    9 101  1  1
#> 51    3    9 102  0  3
#> 52    4    9  97  7  3
#> 53    5    9  88 11 10
#> 54    6    9 102  4  5
#> 55    1   10 100  0  1
#> 56    2   10 100  2  1
#> 57    3   10 102  0  3
#> 58    4   10  97  7  3
#> 59    5   10  87 11 11
#> 60    6   10 101  5  5

## Extract the number of recovered individuals in the first node
## at the time-points in 'tspan'.
trajectory(result, compartments = "R", index = 1)
#>    node time R
#> 1     1    1 0
#> 2     1    2 0
#> 3     1    3 1
#> 4     1    4 1
#> 5     1    5 1
#> 6     1    6 1
#> 7     1    7 1
#> 8     1    8 1
#> 9     1    9 1
#> 10    1   10 1

## Extract the number of recovered individuals in the first and
## third node at the time-points in 'tspan'.
trajectory(result, compartments = "R", index = c(1, 3))
#>    node time R
#> 1     1    1 0
#> 2     3    1 0
#> 3     1    2 0
#> 4     3    2 1
#> 5     1    3 1
#> 6     3    3 1
#> 7     1    4 1
#> 8     3    4 2
#> 9     1    5 1
#> 10    3    5 2
#> 11    1    6 1
#> 12    3    6 3
#> 13    1    7 1
#> 14    3    7 3
#> 15    1    8 1
#> 16    3    8 3
#> 17    1    9 1
#> 18    3    9 3
#> 19    1   10 1
#> 20    3   10 3

## Create an 'SISe' model with 6 nodes and initialize
## it to run over 10 days.
model <- SISe(
  u0 = data.frame(
    S = 100:105,
    I = 1:6),
  tspan = 1:10,
  phi = rep(0, 6),
  upsilon = 0.02,
  gamma = 0.1,
  alpha = 1,
  epsilon = 1.1e-5,
  beta_t1 = 0.15,
  beta_t2 = 0.15,
  beta_t3 = 0.15,
  beta_t4 = 0.15,
  end_t1 = 91,
  end_t2 = 182,
  end_t3 = 273,
  end_t4 = 365
)

## Run the model with a fixed seed for reproducibility.
result <- run(model, seed = 123)

## Extract the continuous state variable 'phi' which represents
## the environmental infectious pressure.
trajectory(result, compartments = "phi")
#>    node time        phi
#> 1     1    1 0.00991199
#> 2     2    1 0.01942848
#> 3     3    1 0.02858243
#> 4     4    1 0.03739418
#> 5     5    1 0.04588256
#> 6     6    1 0.05406505
#> 7     1    2 0.01833718
#> 8     2    2 0.03594268
#> 9     3    2 0.04335368
#> 10    4    2 0.06917923
#> 11    5    2 0.08488274
#> 12    6    2 0.10002035
#> 13    1    3 0.02549859
#> 14    2    3 0.04027102
#> 15    3    3 0.05590925
#> 16    4    3 0.08685073
#> 17    5    3 0.10885857
#> 18    6    3 0.13908235
#> 19    1    4 0.03158580
#> 20    2    4 0.04395010
#> 21    3    4 0.05705767
#> 22    4    4 0.10187150
#> 23    5    4 0.12923803
#> 24    6    4 0.17228505
#> 25    1    5 0.02685893
#> 26    2    5 0.05678606
#> 27    3    5 0.05803383
#> 28    4    5 0.11463916
#> 29    5    5 0.13738627
#> 30    6    5 0.20050735
#> 31    1    6 0.02284109
#> 32    2    6 0.06769663
#> 33    3    6 0.05886357
#> 34    4    6 0.12549167
#> 35    5    6 0.14431226
#> 36    6    6 0.23350531
#> 37    1    7 0.01942592
#> 38    2    7 0.07697061
#> 39    3    7 0.07861646
#> 40    4    7 0.13471630
#> 41    5    7 0.15019936
#> 42    6    7 0.26155358
#> 43    1    8 0.01652304
#> 44    2    8 0.08485349
#> 45    3    8 0.08588261
#> 46    4    8 0.14255724
#> 47    5    8 0.14602908
#> 48    6    8 0.27638559
#> 49    1    9 0.01405558
#> 50    2    9 0.10126268
#> 51    3    9 0.08253503
#> 52    4    9 0.14922204
#> 53    5    9 0.14248434
#> 54    6    9 0.27998380
#> 55    1   10 0.01195824
#> 56    2   10 0.12491923
#> 57    3   10 0.07968958
#> 58    4   10 0.15488711
#> 59    5   10 0.13947131
#> 60    6   10 0.26502426
```
