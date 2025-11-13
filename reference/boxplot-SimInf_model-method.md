# Box plot of number of individuals in each compartment

Produce box-and-whisker plot(s) of the number of individuals in each
model compartment.

## Usage

``` r
# S4 method for class 'SimInf_model'
boxplot(x, compartments = NULL, index = NULL, ...)
```

## Arguments

- x:

  The `model` to plot

- compartments:

  specify the names of the compartments to extract data from. The
  compartments can be specified as a character vector e.g.
  `compartments = c('S', 'I', 'R')`, or as a formula e.g.
  `compartments = ~S+I+R` (see ‘Examples’). Default
  (`compartments=NULL`) includes all compartments.

- index:

  indices specifying the nodes to include when plotting data. Default
  `index = NULL` include all nodes in the model.

- ...:

  Additional arguments affecting the plot produced.

## Examples

``` r
## Create an 'SIR' model with 10 nodes and initialise
## it with 99 susceptible individuals and one infected
## individual. Let the model run over 100 days.
model <- SIR(u0 = data.frame(S = rep(99, 10),
                             I = rep(1, 10),
                             R = rep(0, 10)),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the model and save the result.
result <- run(model)

## Create a boxplot that includes all compartments in all nodes.
boxplot(result)


## Create a boxplot that includes the S and I compartments in
## nodes 1 and 2.
boxplot(result, ~S+I, 1:2)
```
