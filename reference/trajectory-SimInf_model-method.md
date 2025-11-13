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
#> 2     2    1 100  2 1
#> 3     3    1 101  4 0
#> 4     4    1 103  4 0
#> 5     5    1 104  5 0
#> 6     6    1 104  7 0
#> 7     1    2  99  2 0
#> 8     2    2 100  2 1
#> 9     3    2  99  5 1
#> 10    4    2 102  5 0
#> 11    5    2 103  6 0
#> 12    6    2 103  7 1
#> 13    1    3  98  2 1
#> 14    2    3  99  3 1
#> 15    3    3  98  6 1
#> 16    4    3 101  6 0
#> 17    5    3 102  6 1
#> 18    6    3 103  5 3
#> 19    1    4  97  3 1
#> 20    2    4  99  3 1
#> 21    3    4  98  5 2
#> 22    4    4 101  5 1
#> 23    5    4  97  9 3
#> 24    6    4 103  5 3
#> 25    1    5  96  4 1
#> 26    2    5  99  3 1
#> 27    3    5  97  5 3
#> 28    4    5 101  5 1
#> 29    5    5  97  9 3
#> 30    6    5 102  6 3
#> 31    1    6  95  5 1
#> 32    2    6  98  4 1
#> 33    3    6  95  5 5
#> 34    4    6  99  7 1
#> 35    5    6  96 10 3
#> 36    6    6  99  8 4
#> 37    1    7  95  5 1
#> 38    2    7  98  3 2
#> 39    3    7  94  5 6
#> 40    4    7  98  7 2
#> 41    5    7  95  9 5
#> 42    6    7  98  8 5
#> 43    1    8  95  4 2
#> 44    2    8  97  4 2
#> 45    3    8  94  4 7
#> 46    4    8  96  9 2
#> 47    5    8  94  7 8
#> 48    6    8  96 10 5
#> 49    1    9  94  5 2
#> 50    2    9  95  6 2
#> 51    3    9  93  5 7
#> 52    4    9  96  8 3
#> 53    5    9  93  8 8
#> 54    6    9  95  9 7
#> 55    1   10  94  5 2
#> 56    2   10  95  6 2
#> 57    3   10  93  4 8
#> 58    4   10  93 11 3
#> 59    5   10  93  7 9
#> 60    6   10  95  9 7

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
#> 8     1    8 2
#> 9     1    9 2
#> 10    1   10 2

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
#> 10    3    5 3
#> 11    1    6 1
#> 12    3    6 5
#> 13    1    7 1
#> 14    3    7 6
#> 15    1    8 2
#> 16    3    8 7
#> 17    1    9 2
#> 18    3    9 7
#> 19    1   10 2
#> 20    3   10 8

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
#> 5     5    1 0.045882560
#> 6     6    1 0.036047036
#> 7     1    2 0.018337182
#> 8     2    2 0.017981515
#> 9     3    2 0.043353683
#> 10    4    2 0.051889509
#> 11    5    2 0.084882735
#> 12    6    2 0.066687017
#> 13    1    3 0.025498595
#> 14    2    3 0.025004026
#> 15    3    3 0.055909250
#> 16    4    3 0.072154466
#> 17    5    3 0.118032885
#> 18    6    3 0.083721991
#> 19    1    4 0.031585795
#> 20    2    4 0.030973160
#> 21    3    4 0.057057672
#> 22    4    4 0.089379679
#> 23    5    4 0.146210512
#> 24    6    4 0.098201720
#> 25    1    5 0.036759916
#> 26    2    5 0.036046924
#> 27    3    5 0.058033831
#> 28    4    5 0.104021110
#> 29    5    5 0.160987183
#> 30    6    5 0.128527507
#> 31    1    6 0.041157919
#> 32    2    6 0.040359623
#> 33    3    6 0.058863566
#> 34    4    6 0.116466327
#> 35    5    6 0.191895977
#> 36    6    6 0.154304426
#> 37    1    7 0.044896221
#> 38    2    7 0.044025417
#> 39    3    7 0.050045031
#> 40    4    7 0.127044761
#> 41    5    7 0.208994140
#> 42    6    7 0.185223816
#> 43    1    8 0.048073778
#> 44    2    8 0.047141343
#> 45    3    8 0.042549276
#> 46    4    8 0.126690636
#> 47    5    8 0.232701890
#> 48    6    8 0.193487280
#> 49    1    9 0.050774701
#> 50    2    9 0.059498617
#> 51    3    9 0.036177885
#> 52    4    9 0.126389629
#> 53    5    9 0.280376414
#> 54    6    9 0.200511224
#> 55    1   10 0.053070486
#> 56    2   10 0.070002300
#> 57    3   10 0.030762202
#> 58    4   10 0.126133774
#> 59    5   10 0.311725447
#> 60    6   10 0.206481576
```
