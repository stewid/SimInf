# Extract data from a simulated trajectory

Extract the number of individuals in each compartment in every node
after generating a single stochastic trajectory with
[`run`](http://stewid.github.io/SimInf/reference/run.md).

## Usage

``` r
# S4 method for class 'SimInf_model'
trajectory(model, compartments, index, format = c("data.frame", "matrix"))
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

Description of the layout of the internal matrix (`U`) that is returned
if `format = "matrix"`. `U[, j]` contains the number of individuals in
each compartment at `tspan[j]`. `U[1:Nc, j]` contains the number of
individuals in node 1 at `tspan[j]`. `U[(Nc + 1):(2 * Nc), j]` contains
the number of individuals in node 2 at `tspan[j]` etc, where `Nc` is the
number of compartments in the model. The dimension of the matrix is
\\N_n N_c \times\\ `length(tspan)` where \\N_n\\ is the number of nodes.
Since version 10, the internal format of `U` has been expanded to also
allow replicates of each node. This new functionality is used by the
bootstrap filtering algorithm. Each replicate adds new columns to `U` so
that the data for each replicate is in blocks of `length(tspan)`
columns.

## Internal format of the continuous state variables

Description of the layout of the matrix that is returned if
`format = "matrix"`. The result matrix for the real-valued continuous
state. `V[, j]` contains the real-valued state of the system at
`tspan[j]`. The dimension of the matrix is \\N_n\\`dim(ldata)[1]`
\\\times\\ `length(tspan)`. Since version 10, the internal format of `V`
has been expanded to also allow replicates of each node. This new
functionality is used by the bootstrap filtering algorithm. Each
replicate adds new columns to `V` so that the data for each replicate is
in blocks of `length(tspan)` columns.

## Examples

``` r
## Create an 'SIR' model with 6 nodes and initialize
## it to run over 10 days.
u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Extract the number of individuals in each compartment at the
## time-points in 'tspan'.
trajectory(result)
#>    node time   S I R
#> 1     1    1 100 1 0
#> 2     2    1 101 2 0
#> 3     3    1 101 4 0
#> 4     4    1 103 4 0
#> 5     5    1 103 6 0
#> 6     6    1 105 6 0
#> 7     1    2 100 1 0
#> 8     2    2 101 2 0
#> 9     3    2 100 3 2
#> 10    4    2 102 5 0
#> 11    5    2 103 6 0
#> 12    6    2 105 6 0
#> 13    1    3 100 1 0
#> 14    2    3 100 3 0
#> 15    3    3 100 3 2
#> 16    4    3 102 5 0
#> 17    5    3 103 6 0
#> 18    6    3 105 6 0
#> 19    1    4 100 1 0
#> 20    2    4 100 3 0
#> 21    3    4  98 5 2
#> 22    4    4 101 5 1
#> 23    5    4 102 7 0
#> 24    6    4 103 7 1
#> 25    1    5 100 1 0
#> 26    2    5 100 3 0
#> 27    3    5  98 4 3
#> 28    4    5 100 6 1
#> 29    5    5 101 7 1
#> 30    6    5 103 6 2
#> 31    1    6 100 1 0
#> 32    2    6  96 7 0
#> 33    3    6  98 4 3
#> 34    4    6  99 7 1
#> 35    5    6 101 6 2
#> 36    6    6 101 8 2
#> 37    1    7 100 0 1
#> 38    2    7  96 6 1
#> 39    3    7  98 4 3
#> 40    4    7  98 6 3
#> 41    5    7 100 6 3
#> 42    6    7 101 7 3
#> 43    1    8 100 0 1
#> 44    2    8  94 7 2
#> 45    3    8  96 6 3
#> 46    4    8  98 6 3
#> 47    5    8 100 6 3
#> 48    6    8 100 7 4
#> 49    1    9 100 0 1
#> 50    2    9  92 8 3
#> 51    3    9  96 5 4
#> 52    4    9  97 6 4
#> 53    5    9  99 7 3
#> 54    6    9 100 7 4
#> 55    1   10 100 0 1
#> 56    2   10  91 9 3
#> 57    3   10  95 6 4
#> 58    4   10  97 6 4
#> 59    5   10  98 6 5
#> 60    6   10  99 8 4

## Extract the number of recovered individuals in the first node
## at the time-points in 'tspan'.
trajectory(result, compartments = "R", index = 1)
#>    node time R
#> 1     1    1 0
#> 2     1    2 0
#> 3     1    3 0
#> 4     1    4 0
#> 5     1    5 0
#> 6     1    6 0
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
#> 4     3    2 2
#> 5     1    3 0
#> 6     3    3 2
#> 7     1    4 0
#> 8     3    4 2
#> 9     1    5 0
#> 10    3    5 3
#> 11    1    6 0
#> 12    3    6 3
#> 13    1    7 1
#> 14    3    7 3
#> 15    1    8 1
#> 16    3    8 3
#> 17    1    9 1
#> 18    3    9 4
#> 19    1   10 1
#> 20    3   10 4

## Create an 'SISe' model with 6 nodes and initialize
## it to run over 10 days.
u0 <- data.frame(S = 100:105, I = 1:6)
model <- SISe(u0 = u0, tspan = 1:10, phi = rep(0, 6),
    upsilon = 0.02, gamma = 0.1, alpha = 1, epsilon = 1.1e-5,
    beta_t1 = 0.15, beta_t2 = 0.15, beta_t3 = 0.15, beta_t4 = 0.15,
    end_t1 = 91, end_t2 = 182, end_t3 = 273, end_t4 = 365)

## Run the model
result <- run(model)

## Extract the continuous state variable 'phi' which represents
## the environmental infectious pressure.
trajectory(result, "phi")
#>    node time         phi
#> 1     1    1 0.009911990
#> 2     2    1 0.019428476
#> 3     3    1 0.028582429
#> 4     4    1 0.037394178
#> 5     5    1 0.045882560
#> 6     6    1 0.045056045
#> 7     1    2 0.018337182
#> 8     2    2 0.035942680
#> 9     3    2 0.052877493
#> 10    4    2 0.069179229
#> 11    5    2 0.084882735
#> 12    6    2 0.074344674
#> 13    1    3 0.025498595
#> 14    2    3 0.049979754
#> 15    3    3 0.073528298
#> 16    4    3 0.096196522
#> 17    5    3 0.118032885
#> 18    6    3 0.090231000
#> 19    1    4 0.021684805
#> 20    2    4 0.061911266
#> 21    3    4 0.081557672
#> 22    4    4 0.119161221
#> 23    5    4 0.155384824
#> 24    6    4 0.112743386
#> 25    1    5 0.018443085
#> 26    2    5 0.081761790
#> 27    3    5 0.088382640
#> 28    4    5 0.148027010
#> 29    5    5 0.187133972
#> 30    6    5 0.122869905
#> 31    1    6 0.015687622
#> 32    2    6 0.088925997
#> 33    3    6 0.094183863
#> 34    4    6 0.181908725
#> 35    5    6 0.204946435
#> 36    6    6 0.131477447
#> 37    1    7 0.013345479
#> 38    2    7 0.095015573
#> 39    3    7 0.108638712
#> 40    4    7 0.201362388
#> 41    5    7 0.210912718
#> 42    6    7 0.138793857
#> 43    1    8 0.011354657
#> 44    2    8 0.090482975
#> 45    3    8 0.111401524
#> 46    4    8 0.227243796
#> 47    5    8 0.215984058
#> 48    6    8 0.154021814
#> 49    1    9 0.009662458
#> 50    2    9 0.086630267
#> 51    3    9 0.113749915
#> 52    4    9 0.249242993
#> 53    5    9 0.220294697
#> 54    6    9 0.166965578
#> 55    1   10 0.008224090
#> 56    2   10 0.083355465
#> 57    3   10 0.115746047
#> 58    4   10 0.267942310
#> 59    5   10 0.223958740
#> 60    6   10 0.186976786
```
