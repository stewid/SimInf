# Example data to initialize the ‘SEIR’ model

Example data to initialize a population of 1600 nodes and demonstrate
the [`SEIR`](http://stewid.github.io/SimInf/reference/SEIR-class.md)
model.

## Usage

``` r
u0_SEIR()
```

## Value

A `data.frame`

## Details

A `data.frame` with the number of individuals in the ‘S’, ‘E’, ‘I’ and
‘R’ compartments in 1600 nodes. Note that the ‘E’, ‘I’ and ‘R’
compartments are zero.

## Examples

``` r
if (FALSE) { # \dontrun{
## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SEIR' model with 1600 nodes and initialize it to
## run over 4*365 days and record data at weekly time-points.
## Add ten infected individuals to the first node.
u0 <- u0_SEIR()
u0$I[1] <- 10
tspan <- seq(from = 1, to = 4*365, by = 7)
model <- SEIR(u0      = u0,
              tspan   = tspan,
              events  = events_SEIR(),
              beta    = 0.16,
              epsilon = 0.25,
              gamma   = 0.01)

## Run the model to generate a single stochastic trajectory.
result <- run(model)
plot(result)

## Summarize trajectory
summary(result)
} # }
```
