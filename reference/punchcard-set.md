# Set a template for where to record result during a simulation

Using a sparse result matrix can save a lot of memory if the model
contains many nodes and time-points, but where only a few of the data
points are of interest for post-processing.

## Usage

``` r
punchcard(model) <- value

# S4 method for class 'SimInf_model'
punchcard(model) <- value
```

## Arguments

- model:

  The `model` to set a template for where to record result.

- value:

  A `data.frame` that specify the nodes, time-points and compartments to
  record the number of individuals at `tspan`. Use `NULL` to reset the
  model to record the number of inidividuals in each compartment in
  every node at each time-point in tspan.

## Details

Using a sparse result matrix can save a lot of memory if the model
contains many nodes and time-points, but where only a few of the data
points are of interest for post-processing. To use this feature, a
template has to be defined for which data points to record. This is done
using a `data.frame` that specifies the time-points (column ‘time’) and
nodes (column ‘node’) to record the state of the compartments, see
‘Examples’. The specified time-points, nodes and compartments must exist
in the model, or an error is raised. Note that specifying a template
only affects which data-points are recorded for post-processing, it does
not affect how the solver simulates the trajectory.

## Examples

``` r
## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 6 nodes and initialize it to run over 10 days.
u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Run the model.
result <- run(model)

## Display the trajectory with data for every node at each
## time-point in tspan.
trajectory(result)
#>    node time   S  I R
#> 1     1    1 100  1 0
#> 2     2    1 101  2 0
#> 3     3    1 102  3 0
#> 4     4    1 102  5 0
#> 5     5    1 103  6 0
#> 6     6    1 105  6 0
#> 7     1    2 100  1 0
#> 8     2    2 101  2 0
#> 9     3    2 101  4 0
#> 10    4    2 101  5 1
#> 11    5    2 103  6 0
#> 12    6    2 105  6 0
#> 13    1    3  99  2 0
#> 14    2    3 101  2 0
#> 15    3    3 101  4 0
#> 16    4    3  99  6 2
#> 17    5    3 101  8 0
#> 18    6    3 103  7 1
#> 19    1    4  98  3 0
#> 20    2    4 101  2 0
#> 21    3    4 101  4 0
#> 22    4    4  98  6 3
#> 23    5    4  99 10 0
#> 24    6    4 101  8 2
#> 25    1    5  98  3 0
#> 26    2    5 101  2 0
#> 27    3    5 100  5 0
#> 28    4    5  97  6 4
#> 29    5    5  98  9 2
#> 30    6    5 101  6 4
#> 31    1    6  98  2 1
#> 32    2    6 101  2 0
#> 33    3    6 100  5 0
#> 34    4    6  97  5 5
#> 35    5    6  98  8 3
#> 36    6    6 100  7 4
#> 37    1    7  98  2 1
#> 38    2    7  98  5 0
#> 39    3    7 100  5 0
#> 40    4    7  92 10 5
#> 41    5    7  98  7 4
#> 42    6    7  99  8 4
#> 43    1    8  97  3 1
#> 44    2    8  98  5 0
#> 45    3    8  98  6 1
#> 46    4    8  92  8 7
#> 47    5    8  95 10 4
#> 48    6    8  99  8 4
#> 49    1    9  97  3 1
#> 50    2    9  97  6 0
#> 51    3    9  98  4 3
#> 52    4    9  91  9 7
#> 53    5    9  94 10 5
#> 54    6    9  99  7 5
#> 55    1   10  97  3 1
#> 56    2   10  96  6 1
#> 57    3   10  98  4 3
#> 58    4   10  89 11 7
#> 59    5   10  93  9 7
#> 60    6   10  98  8 5

## Assume we are only interested in nodes '2' and '4' at the
## time-points '3' and '5'
df <- data.frame(time = c(3, 5, 3, 5),
                 node = c(2, 2, 4, 4),
                 S = c(TRUE, TRUE, TRUE, TRUE),
                 I = c(TRUE, TRUE, TRUE, TRUE),
                 R = c(TRUE, TRUE, TRUE, TRUE))
punchcard(model) <- df
result <- run(model)
trajectory(result)
#>   node time   S I R
#> 1    2    3 100 3 0
#> 2    4    3 102 5 0
#> 3    2    5 100 3 0
#> 4    4    5 100 6 1

## We can also specify to record only some of the compartments in
## each time-step.
df <- data.frame(time = c(3, 5, 3, 5),
                 node = c(2, 2, 4, 4),
                 S = c(FALSE, TRUE, TRUE, TRUE),
                 I = c(TRUE, FALSE, TRUE, FALSE),
                 R = c(TRUE, FALSE, TRUE, TRUE))
punchcard(model) <- df
result <- run(model)
trajectory(result)
#>   node time  S  I  R
#> 1    2    3 NA  4  0
#> 2    4    3 99  7  1
#> 3    2    5 98 NA NA
#> 4    4    5 99 NA  3

## A shortcut to specify to record all of the compartments in
## each time-step is to only inlude node and time.
df <- data.frame(time = c(3, 5, 3, 5),
                 node = c(2, 2, 4, 4))
punchcard(model) <- df
result <- run(model)
trajectory(result)
#>   node time   S I R
#> 1    2    3  97 6 0
#> 2    4    3 103 4 0
#> 3    2    5  96 6 1
#> 4    4    5 101 6 0

## It is possible to use an empty 'data.frame' to specify
## that no data-points should be recorded for the trajectory.
punchcard(model) <- data.frame()
result <- run(model)
trajectory(result)
#> [1] node time S    I    R   
#> <0 rows> (or 0-length row.names)

## Use 'NULL' to reset the model to record data for every node at
## each time-point in tspan.
punchcard(model) <- NULL
result <- run(model)
trajectory(result)
#>    node time   S  I R
#> 1     1    1 100  1 0
#> 2     2    1 101  2 0
#> 3     3    1  99  6 0
#> 4     4    1 103  4 0
#> 5     5    1 104  5 0
#> 6     6    1 105  6 0
#> 7     1    2 100  1 0
#> 8     2    2 100  3 0
#> 9     3    2  98  7 0
#> 10    4    2 102  4 1
#> 11    5    2 104  5 0
#> 12    6    2 104  7 0
#> 13    1    3 100  1 0
#> 14    2    3 100  2 1
#> 15    3    3  98  7 0
#> 16    4    3 101  5 1
#> 17    5    3 101  8 0
#> 18    6    3 103  6 2
#> 19    1    4 100  1 0
#> 20    2    4 100  2 1
#> 21    3    4  96  9 0
#> 22    4    4  99  7 1
#> 23    5    4 101  6 2
#> 24    6    4 103  5 3
#> 25    1    5 100  1 0
#> 26    2    5 100  2 1
#> 27    3    5  96  7 2
#> 28    4    5  99  6 2
#> 29    5    5  96 11 2
#> 30    6    5 103  5 3
#> 31    1    6 100  1 0
#> 32    2    6  99  3 1
#> 33    3    6  95  8 2
#> 34    4    6  98  7 2
#> 35    5    6  95 11 3
#> 36    6    6 102  5 4
#> 37    1    7 100  1 0
#> 38    2    7  99  3 1
#> 39    3    7  94  9 2
#> 40    4    7  98  7 2
#> 41    5    7  94 12 3
#> 42    6    7 102  4 5
#> 43    1    8 100  1 0
#> 44    2    8  99  3 1
#> 45    3    8  94  8 3
#> 46    4    8  97  7 3
#> 47    5    8  93 13 3
#> 48    6    8 101  5 5
#> 49    1    9 100  1 0
#> 50    2    9  99  3 1
#> 51    3    9  90 11 4
#> 52    4    9  97  7 3
#> 53    5    9  91 13 5
#> 54    6    9 101  5 5
#> 55    1   10 100  1 0
#> 56    2   10  99  3 1
#> 57    3   10  87 14 4
#> 58    4   10  97  5 5
#> 59    5   10  89 14 6
#> 60    6   10  99  7 5
```
