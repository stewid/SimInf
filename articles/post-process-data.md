# Post-process data in a trajectory

After a model is created, a simulation is started with a call to the
[`run()`](http://stewid.github.io/SimInf/reference/run.md) function. The
function returns a modified model object containing a single stochastic
solution trajectory. This trajectory includes the state of each
compartment recorded at every time-point specified in `tspan`.

This vignette introduces the functionality in `SimInf` to post-process
and explore this trajectory data.

Let us first load the SimInf package.

``` r

library(SimInf)
```

## Extract trajectory data with `trajectory()`

Most modeling studies require custom data analysis beyond simple
plotting. To support this, `SimInf` provides the
[`trajectory()`](http://stewid.github.io/SimInf/reference/trajectory.md)
method to extract the raw data as a `data.frame`. This is useful if you
need to:

- Perform custom statistical calculations (e.g., time to peak).
- Export data to CSV for use in other software.
- Combine results from multiple simulation runs.

Let’s simulate 10 days of data from an SIR model with 6 nodes. For
reproducibility, we set the number of threads and specify a seed.

``` r

set_num_threads(1)

u0 <- data.frame(
  S = c(100, 101, 102, 103, 104, 105),
  I = c(1, 2, 3, 4, 5, 6),
  R = c(0, 0, 0, 0, 0, 0)
)

model  <- SIR(
  u0 = u0,
  tspan = 1:10,
  beta = 0.16,
  gamma = 0.077
)

result <- run(model, seed = 123)
```

Extract the full trajectory data (all compartments, all nodes).

``` r

trajectory(result)
```

    ##    node time   S  I  R
    ## 1     1    1 100  1  0
    ## 2     2    1 101  2  0
    ## 3     3    1 102  3  0
    ## 4     4    1 103  4  0
    ## 5     5    1 100  9  0
    ## 6     6    1 104  7  0
    ## 7     1    2 100  1  0
    ## 8     2    2 101  2  0
    ## 9     3    2 102  2  1
    ## 10    4    2 103  3  1
    ## 11    5    2  99  9  1
    ## 12    6    2 104  6  1
    ## 13    1    3 100  0  1
    ## 14    2    3 101  2  0
    ## 15    3    3 102  2  1
    ## 16    4    3 103  3  1
    ## 17    5    3  97  8  4
    ## 18    6    3 104  6  1
    ## 19    1    4 100  0  1
    ## 20    2    4 101  2  0
    ## 21    3    4 102  1  2
    ## 22    4    4 103  3  1
    ## 23    5    4  95 10  4
    ## 24    6    4 103  7  1
    ## 25    1    5 100  0  1
    ## 26    2    5 101  2  0
    ## 27    3    5 102  1  2
    ## 28    4    5 101  5  1
    ## 29    5    5  94 11  4
    ## 30    6    5 103  7  1
    ## 31    1    6 100  0  1
    ## 32    2    6 101  2  0
    ## 33    3    6 102  0  3
    ## 34    4    6 100  6  1
    ## 35    5    6  92 11  6
    ## 36    6    6 102  7  2
    ## 37    1    7 100  0  1
    ## 38    2    7 101  1  1
    ## 39    3    7 102  0  3
    ## 40    4    7 100  5  2
    ## 41    5    7  89 13  7
    ## 42    6    7 102  5  4
    ## 43    1    8 100  0  1
    ## 44    2    8 101  1  1
    ## 45    3    8 102  0  3
    ## 46    4    8  98  7  2
    ## 47    5    8  88 14  7
    ## 48    6    8 102  5  4
    ## 49    1    9 100  0  1
    ## 50    2    9 101  1  1
    ## 51    3    9 102  0  3
    ## 52    4    9  97  7  3
    ## 53    5    9  88 11 10
    ## 54    6    9 102  4  5
    ## 55    1   10 100  0  1
    ## 56    2   10 100  2  1
    ## 57    3   10 102  0  3
    ## 58    4   10  97  7  3
    ## 59    5   10  87 11 11
    ## 60    6   10 101  5  5

Extract the number of recovered individuals (R) in the first node only.

``` r

trajectory(result, compartments = "R", index = 1)
```

    ##    node time R
    ## 1     1    1 0
    ## 2     1    2 0
    ## 3     1    3 1
    ## 4     1    4 1
    ## 5     1    5 1
    ## 6     1    6 1
    ## 7     1    7 1
    ## 8     1    8 1
    ## 9     1    9 1
    ## 10    1   10 1

Extract the number of recovered individuals in the first and third
nodes.

``` r

trajectory(result, compartments = "R", index = c(1, 3))
```

    ##    node time R
    ## 1     1    1 0
    ## 2     3    1 0
    ## 3     1    2 0
    ## 4     3    2 1
    ## 5     1    3 1
    ## 6     3    3 1
    ## 7     1    4 1
    ## 8     3    4 2
    ## 9     1    5 1
    ## 10    3    5 2
    ## 11    1    6 1
    ## 12    3    6 3
    ## 13    1    7 1
    ## 14    3    7 3
    ## 15    1    8 1
    ## 16    3    8 3
    ## 17    1    9 1
    ## 18    3    9 3
    ## 19    1   10 1
    ## 20    3   10 3

## Calculate prevalence from a trajectory using `prevalence()`

The
[`prevalence()`](http://stewid.github.io/SimInf/reference/prevalence.md)
function calculates the proportion of individuals with the disease. It
takes a model object and a formula:

- **Left-hand side (LHS):** Compartments representing “cases” (e.g., I).
- **Right-hand side (RHS):** Compartments representing the “at-risk”
  population (e.g., S + I + R).

The function also supports a `level` argument to change the aggregation
level:

- `level = 1` (default): Prevalence aggregated over all nodes (global).
- `level = 2`: Proportion of nodes that have at least one case.
- `level = 3`: Prevalence calculated within each node (returns a
  matrix).

Let’s determine the proportion of infected individuals in the total
population.

``` r

prevalence(result, I ~ S + I + R)
```

    ##    time prevalence
    ## 1     1 0.04088050
    ## 2     2 0.03616352
    ## 3     3 0.03301887
    ## 4     4 0.03616352
    ## 5     5 0.04088050
    ## 6     6 0.04088050
    ## 7     7 0.03773585
    ## 8     8 0.04245283
    ## 9     9 0.03616352
    ## 10   10 0.03930818

Identical result is obtained with the shorthand `I ~ .` (where `.` means
“all compartments”).

``` r

prevalence(result, I ~ .)
```

    ##    time prevalence
    ## 1     1 0.04088050
    ## 2     2 0.03616352
    ## 3     3 0.03301887
    ## 4     4 0.03616352
    ## 5     5 0.04088050
    ## 6     6 0.04088050
    ## 7     7 0.03773585
    ## 8     8 0.04245283
    ## 9     9 0.03616352
    ## 10   10 0.03930818

Calculate the proportion of nodes that are infected (at least one I
individual).

``` r

prevalence(result, I ~ S + I + R, level = 2)
```

    ##    time prevalence
    ## 1     1  1.0000000
    ## 2     2  1.0000000
    ## 3     3  0.8333333
    ## 4     4  0.8333333
    ## 5     5  0.8333333
    ## 6     6  0.6666667
    ## 7     7  0.6666667
    ## 8     8  0.6666667
    ## 9     9  0.6666667
    ## 10   10  0.6666667

Calculate the prevalence **within each node** individually.

``` r

prevalence(result, I ~ S + I + R, level = 3)
```

    ##    node time  prevalence
    ## 1     1    1 0.009900990
    ## 2     2    1 0.019417476
    ## 3     3    1 0.028571429
    ## 4     4    1 0.037383178
    ## 5     5    1 0.082568807
    ## 6     6    1 0.063063063
    ## 7     1    2 0.009900990
    ## 8     2    2 0.019417476
    ## 9     3    2 0.019047619
    ## 10    4    2 0.028037383
    ## 11    5    2 0.082568807
    ## 12    6    2 0.054054054
    ## 13    1    3 0.000000000
    ## 14    2    3 0.019417476
    ## 15    3    3 0.019047619
    ## 16    4    3 0.028037383
    ## 17    5    3 0.073394495
    ## 18    6    3 0.054054054
    ## 19    1    4 0.000000000
    ## 20    2    4 0.019417476
    ## 21    3    4 0.009523810
    ## 22    4    4 0.028037383
    ## 23    5    4 0.091743119
    ## 24    6    4 0.063063063
    ## 25    1    5 0.000000000
    ## 26    2    5 0.019417476
    ## 27    3    5 0.009523810
    ## 28    4    5 0.046728972
    ## 29    5    5 0.100917431
    ## 30    6    5 0.063063063
    ## 31    1    6 0.000000000
    ## 32    2    6 0.019417476
    ## 33    3    6 0.000000000
    ## 34    4    6 0.056074766
    ## 35    5    6 0.100917431
    ## 36    6    6 0.063063063
    ## 37    1    7 0.000000000
    ## 38    2    7 0.009708738
    ## 39    3    7 0.000000000
    ## 40    4    7 0.046728972
    ## 41    5    7 0.119266055
    ## 42    6    7 0.045045045
    ## 43    1    8 0.000000000
    ## 44    2    8 0.009708738
    ## 45    3    8 0.000000000
    ## 46    4    8 0.065420561
    ## 47    5    8 0.128440367
    ## 48    6    8 0.045045045
    ## 49    1    9 0.000000000
    ## 50    2    9 0.009708738
    ## 51    3    9 0.000000000
    ## 52    4    9 0.065420561
    ## 53    5    9 0.100917431
    ## 54    6    9 0.036036036
    ## 55    1   10 0.000000000
    ## 56    2   10 0.019417476
    ## 57    3   10 0.000000000
    ## 58    4   10 0.065420561
    ## 59    5   10 0.100917431
    ## 60    6   10 0.045045045

## Visualize a trajectory with `plot()`

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) function
provides a quick way to inspect the outcome. It can display:

- The median and quantile range across all nodes.
- Individual trajectories for specific nodes.
- Prevalence curves.

*Note: Since the simulation is stochastic, the exact lines shown below
will vary unless set.seed() is used.*

### Aggregated View (Median and Range)

Plot the median and interquartile range (IQR) of all compartments.

``` r

plot(result)
```

![](post-process-data_files/figure-html/unnamed-chunk-10-1.png)

Plot the median and the middle 95% quantile range.

``` r

plot(result, range = 0.95)
```

![](post-process-data_files/figure-html/unnamed-chunk-11-1.png)

Plot only the infected individuals (I).

``` r

plot(result, "I")
```

![](post-process-data_files/figure-html/unnamed-chunk-12-1.png)

Use formula notation to plot the infected individuals.

``` r

plot(result, ~I)
```

![](post-process-data_files/figure-html/unnamed-chunk-13-1.png)

### Individual Node View

Plot the trajectories for the first three nodes. We use `range = FALSE`
to suppress the shaded median/range bands and show the individual lines.

``` r

plot(result, index = 1:3, range = FALSE)
```

![](post-process-data_files/figure-html/unnamed-chunk-14-1.png)

Use `type = "l"` to draw a line.

``` r

plot(result, index = 1:3, range = FALSE, type = "l")
```

![](post-process-data_files/figure-html/unnamed-chunk-15-1.png)

Plot the infected individuals in the first node only.

``` r

plot(result, "I", index = 1, range = FALSE)
```

![](post-process-data_files/figure-html/unnamed-chunk-16-1.png)

### Prevalence Plots

Plot the proportion of infected individuals in the population.

``` r

plot(result, I ~ S + I + R)
```

![](post-process-data_files/figure-html/unnamed-chunk-17-1.png)

Plot the proportion of nodes with infected individuals (`level = 2`).

``` r

plot(result, I ~ S + I + R, level = 2)
```

![](post-process-data_files/figure-html/unnamed-chunk-18-1.png)

Plot the median and IQR of the prevalence **within in each node**
(`level = 3`).

``` r

plot(result, I ~ S + I + R, level = 3)
```

![](post-process-data_files/figure-html/unnamed-chunk-19-1.png)

Plot the prevalence in the first three nodes.

``` r

plot(result, I ~ S + I + R, level = 3, index = 1:3, range = FALSE)
```

![](post-process-data_files/figure-html/unnamed-chunk-20-1.png)

## Summary

- Use
  [`trajectory()`](http://stewid.github.io/SimInf/reference/trajectory.md)
  to extract raw data for custom analysis.
- Use
  [`prevalence()`](http://stewid.github.io/SimInf/reference/prevalence.md)
  to calculate disease proportions at different aggregation levels.
- Use [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for quick
  visual inspection of medians, ranges, or individual trajectories.

To find more details on the plot method for `SimInf_model` objects, run:

``` r

help("plot,SimInf_model-method", package = "SimInf")
```
