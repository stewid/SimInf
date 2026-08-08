# Brief summary of a `SimInf_model` object

Display key characteristics of a
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object, including the model name, number of nodes, number of replicates
(if greater than one), number of transitions, scheduled events, global
data, local data, continuous state variables, and compartments.

## Usage

``` r
# S4 method for class 'SimInf_model'
show(object)
```

## Arguments

- object:

  The
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  object to display.

## Value

The `object`, returned invisibly.

## Details

The output differs depending on whether the model has been run: before
running, the result matrices (`U` and `V`) are empty; after calling
[`run`](http://stewid.github.io/SimInf/reference/run.md), they contain
the simulated trajectory data and the summary reflects the computed
results.

## See also

[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model.md)
for creating model objects,
[`run`](http://stewid.github.io/SimInf/reference/run.md) for simulating
a trajectory from a model,
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for the class definition.

## Examples

``` r
## For reproducibility, set the seed and number of threads.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 10 nodes and initialise
## it to run over 100 days.
model <- SIR(
  u0 = data.frame(
    S = rep(99, 10),
    I = rep(1, 10),
    R = rep(0, 10)
  ),
  tspan = 1:100,
  beta = 0.16,
  gamma = 0.077
)

## Brief summary of the model
model
#> Model: SIR
#> Number of nodes: 10
#> Number of transitions: 2
#> Number of scheduled events: 0
#> 
#> Local data
#> ----------
#>  Parameter Value
#>  beta      0.160
#>  gamma     0.077
#> 
#> Compartments
#> ------------
#>  - Empty, please run the model first

## Run the model and save the result
result <- run(model)

## Brief summary of the result. Note that 'U' and 'V' are
## non-empty after running the model.
result
#> Model: SIR
#> Number of nodes: 10
#> Number of transitions: 2
#> Number of scheduled events: 0
#> 
#> Local data
#> ----------
#>  Parameter Value
#>  beta      0.160
#>  gamma     0.077
#> 
#> Compartments
#> ------------
#>     Min. 1st Qu. Median  Mean 3rd Qu.  Max.
#>  S  8.00   30.00  81.50 66.17   99.00 99.00
#>  I  0.00    0.00   4.50  7.67   15.00 33.00
#>  R  0.00    1.00   8.00 26.17   57.00 87.00
```
