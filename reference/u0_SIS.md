# Example Initial Population Data for the SIS Model

Dataset containing the initial number of susceptible and infected cattle
across 1,600 herds. Provides realistic population structure for
demonstrating SIS model simulations in a cattle disease epidemiology
context.

## Usage

``` r
u0_SIS()
```

## Value

A `data.frame` with 1,600 rows (one per herd) and 2 columns:

- S:

  Number of susceptible cattle in the herd

- I:

  Number of infected cattle in the herd (all zero at start)

## Details

This dataset represents initial disease states in a population of 1,600
cattle herds (nodes). Each row represents a single herd (node), derived
from the cattle population data by extracting susceptible and infected
compartments. The SIS model is appropriate for diseases where recovered
individuals do not gain immunity.

The data contains:

- S:

  Total susceptible cattle in the herd

- I:

  Total infected cattle (initialized to zero)

The herd size distribution reflects realistic heterogeneity observed in
cattle populations, making it suitable for testing disease dynamics in
endemic situations where immunity is temporary or absent.

## See also

[`SIS`](http://stewid.github.io/SimInf/reference/SIS.md) for creating
SIS models with this initial state and
[`events_SIS`](http://stewid.github.io/SimInf/reference/events_SIS.md)
for associated cattle movement and demographic events

## Examples

``` r
if (FALSE) { # \dontrun{
## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIS' model with 1600 cattle herds (nodes) and
## initialize it to run over 4*365 days. Add one infected animal
## to the first herd to seed the outbreak. Define 'tspan' to record
## the state of the system at daily time-points. Load scheduled
## events for the population of nodes with births, deaths and
## between-node movements of individuals.
u0 <- u0_SIS()
u0$I[1] <- 1
model <- SIS(u0     = u0,
             tspan  = seq(from = 1, to = 4*365, by = 1),
             events = events_SIS(),
             beta   = 0.16,
             gamma  = 0.01)

## Run the model to generate a single stochastic trajectory.
result <- run(model)
plot(result)

## Plot the trajectory for the first herd.
plot(result, index = 1)

## Summarize trajectory
summary(result)
} # }
```
