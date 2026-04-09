# Box plot of number of individuals in each compartment

Produce box-and-whisker plots summarizing the distribution of the number
of individuals in specified model compartments. The plots aggregate data
across all time points in `tspan` and the selected nodes (specified by
`index`).

## Usage

``` r
# S4 method for class 'SimInf_model'
boxplot(x, compartments = NULL, index = NULL, ...)
```

## Arguments

- x:

  The `SimInf_model` object containing the trajectory data.

- compartments:

  Specify the names of the compartments to include. Can be a character
  vector (e.g., `c("S", "I", "R")`) or a formula (e.g., `~S + I + R`).
  If `NULL` (default), all compartments in the model are included.

- index:

  Indices of the nodes to include in the plot. If `NULL` (default), data
  from all nodes are pooled together. If a vector (e.g., `1:2`), only
  data from those specific nodes are used.

- ...:

  Additional graphical arguments passed to
  [`boxplot()`](https://rdrr.io/r/graphics/boxplot.html).

## Details

This function is useful for visualizing the variability and range of
compartment counts across the entire simulation trajectory.

## Examples

``` r
## For reproducibility, set the seed and number of threads.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 10 nodes.
model <- SIR(u0 = data.frame(S = rep(99, 10),
                             I = rep(1, 10),
                             R = rep(0, 10)),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the model.
result <- run(model)

## Create a boxplot for all compartments across all nodes and time
## points.
boxplot(result)


## Create a boxplot for the S and I compartments, but only for
## nodes 1 and 2.
boxplot(result, compartments = c("S", "I"), index = 1:2)
```
