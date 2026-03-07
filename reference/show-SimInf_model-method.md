# Brief summary of `SimInf_model`

Brief summary of `SimInf_model`

## Usage

``` r
# S4 method for class 'SimInf_model'
show(object)
```

## Arguments

- object:

  The SimInf_model `object`

## Value

None (invisible 'NULL').

## Examples

``` r
## Create an 'SIR' model with 10 nodes and initialise
## it to run over 100 days.
model <- SIR(u0 = data.frame(S = rep(99, 10),
                             I = rep(1, 10),
                             R = rep(0, 10)),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

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
#>  S 18.00   68.00  99.00 81.09   99.00 99.00
#>  I  0.00    0.00   0.00  3.97    6.00 26.00
#>  R  0.00    1.00   1.00 14.94   18.00 78.00
```
