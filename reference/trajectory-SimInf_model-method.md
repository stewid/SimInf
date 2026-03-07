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
#>    node time   S  I R
#> 1     1    1 100  1 0
#> 2     2    1 100  3 0
#> 3     3    1 102  3 0
#> 4     4    1 102  5 0
#> 5     5    1 102  7 0
#> 6     6    1 103  8 0
#> 7     1    2 100  1 0
#> 8     2    2 100  3 0
#> 9     3    2 102  3 0
#> 10    4    2 101  5 1
#> 11    5    2 100  8 1
#> 12    6    2 103  8 0
#> 13    1    3 100  1 0
#> 14    2    3  99  4 0
#> 15    3    3 102  2 1
#> 16    4    3  99  7 1
#> 17    5    3 100  7 2
#> 18    6    3 101 10 0
#> 19    1    4 100  1 0
#> 20    2    4  98  5 0
#> 21    3    4 102  1 2
#> 22    4    4  99  5 3
#> 23    5    4 100  5 4
#> 24    6    4 100  9 2
#> 25    1    5 100  1 0
#> 26    2    5  98  4 1
#> 27    3    5 102  0 3
#> 28    4    5  99  5 3
#> 29    5    5  99  6 4
#> 30    6    5 100  8 3
#> 31    1    6 100  1 0
#> 32    2    6  97  5 1
#> 33    3    6 102  0 3
#> 34    4    6  99  5 3
#> 35    5    6  99  6 4
#> 36    6    6  99  9 3
#> 37    1    7 100  1 0
#> 38    2    7  97  5 1
#> 39    3    7 102  0 3
#> 40    4    7  95  8 4
#> 41    5    7  99  6 4
#> 42    6    7  96 12 3
#> 43    1    8 100  1 0
#> 44    2    8  96  5 2
#> 45    3    8 102  0 3
#> 46    4    8  95  8 4
#> 47    5    8  98  7 4
#> 48    6    8  95 11 5
#> 49    1    9 100  1 0
#> 50    2    9  95  6 2
#> 51    3    9 102  0 3
#> 52    4    9  94  7 6
#> 53    5    9  97  8 4
#> 54    6    9  91 15 5
#> 55    1   10 100  1 0
#> 56    2   10  94  7 2
#> 57    3   10 102  0 3
#> 58    4   10  94  6 7
#> 59    5   10  95 10 4
#> 60    6   10  88 18 5

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
#> 7     1    7 0
#> 8     1    8 0
#> 9     1    9 0
#> 10    1   10 0

## Extract the number of recovered individuals in the first and
## third node at the time-points in 'tspan'.
trajectory(result, compartments = "R", index = c(1, 3))
#>    node time R
#> 1     1    1 0
#> 2     3    1 0
#> 3     1    2 0
#> 4     3    2 0
#> 5     1    3 0
#> 6     3    3 1
#> 7     1    4 0
#> 8     3    4 2
#> 9     1    5 0
#> 10    3    5 3
#> 11    1    6 0
#> 12    3    6 3
#> 13    1    7 0
#> 14    3    7 3
#> 15    1    8 0
#> 16    3    8 3
#> 17    1    9 0
#> 18    3    9 3
#> 19    1   10 0
#> 20    3   10 3

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
#> 2     2    1 0.009719738
#> 3     3    1 0.028582429
#> 4     4    1 0.028048383
#> 5     5    1 0.036708248
#> 6     6    1 0.054065054
#> 7     1    2 0.008436192
#> 8     2    2 0.017981515
#> 9     3    2 0.052877493
#> 10    4    2 0.051889509
#> 11    5    2 0.067910258
#> 12    6    2 0.100020350
#> 13    1    3 0.007181763
#> 14    2    3 0.025004026
#> 15    3    3 0.064004488
#> 16    4    3 0.062808671
#> 17    5    3 0.094431967
#> 18    6    3 0.130073343
#> 19    1    4 0.006115498
#> 20    2    4 0.030973160
#> 21    3    4 0.073462434
#> 22    4    4 0.072089959
#> 23    5    4 0.116975420
#> 24    6    4 0.155618386
#> 25    1    5 0.005209174
#> 26    2    5 0.036046924
#> 27    3    5 0.081501688
#> 28    4    5 0.089324849
#> 29    5    5 0.126963043
#> 30    6    5 0.177331673
#> 31    1    6 0.004438798
#> 32    2    6 0.040359623
#> 33    3    6 0.088335054
#> 34    4    6 0.094628710
#> 35    5    6 0.144626834
#> 36    6    6 0.195787967
#> 37    1    7 0.003783978
#> 38    2    7 0.044025417
#> 39    3    7 0.094143415
#> 40    4    7 0.108482787
#> 41    5    7 0.159641057
#> 42    6    7 0.220484826
#> 43    1    8 0.003227381
#> 44    2    8 0.047141343
#> 45    3    8 0.099080522
#> 46    4    8 0.120258752
#> 47    5    8 0.172403146
#> 48    6    8 0.241477156
#> 49    1    9 0.002754274
#> 50    2    9 0.049789879
#> 51    3    9 0.103277062
#> 52    4    9 0.130268322
#> 53    5    9 0.174076610
#> 54    6    9 0.277338655
#> 55    1   10 0.002352133
#> 56    2   10 0.052041135
#> 57    3   10 0.106844122
#> 58    4   10 0.138776457
#> 59    5   10 0.184673366
#> 60    6   10 0.298811920
```
