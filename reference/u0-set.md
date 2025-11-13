# Update the initial compartment state u0 in each node

Update the initial compartment state u0 in each node

## Usage

``` r
u0(model) <- value

# S4 method for class 'SimInf_model'
u0(model) <- value
```

## Arguments

- model:

  The model to update the initial compartment state `u0`.

- value:

  A `data.frame` with the initial state in each node. Each row is one
  node, and the number of rows in `u0` must match the number of nodes in
  `model`. Only the columns in `u0` with a name that matches a
  compartment in the `model` will be used.

## Examples

``` r
## Create an SIR model object.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the SIR model and plot the result.
set.seed(22)
result <- run(model)
plot(result)


## Update u0 and run the model again
u0(model) <- data.frame(S = 990, I = 10, R = 0)
result <- run(model)
plot(result)
```
