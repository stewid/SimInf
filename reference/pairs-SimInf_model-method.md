# Scatterplot matrix of number of individuals in each compartment

Produce a matrix of scatterplots showing the relationship between the
number of individuals in different model compartments. The `ij`th panel
displays the counts of compartment `j` plotted against compartment `i`.

## Usage

``` r
# S4 method for class 'SimInf_model'
pairs(x, compartments = NULL, index = NULL, ...)
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
  [`pairs()`](https://rdrr.io/r/graphics/pairs.html).

## Details

Data is aggregated across all time points in `tspan` and the selected
nodes (specified by `index`). This allows for visualizing the
correlation structure between compartments throughout the simulation.

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

## Create a scatterplot matrix for all compartments across all
## nodes.
pairs(result)


## Create a scatterplot matrix for the S and I compartments,
## using only data from nodes 1 and 2.
pairs(result, compartments = c("S", "I"), index = 1:2)
```
