# Example data to initialize the ‘SIS’ model

Example data to initialize a population of 1600 nodes and demonstrate
the [`SIS`](http://stewid.github.io/SimInf/reference/SIS-class.md)
model.

## Usage

``` r
u0_SIS()
```

## Value

A `data.frame`

## Details

A `data.frame` with the number of individuals in the ‘S’, and ‘I’
compartments in 1600 nodes. Note that the ‘I’ compartment is zero.

## Examples

``` r
if (FALSE) { # \dontrun{
## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIS' model with 1600 nodes and initialize
## it to run over 4*365 days. Add one infected individual
## to the first node.
u0 <- u0_SIS()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 4*365, by = 1)
model <- SIS(u0     = u0,
             tspan  = tspan,
             events = events_SIS(),
             beta   = 0.16,
             gamma  = 0.01)

## Run the model to generate a single stochastic trajectory.
result <- run(model)
plot(result)

## Summarize trajectory
summary(result)
} # }
```
