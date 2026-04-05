# Example initial population data for the SISe model

Dataset containing the initial number of susceptible and infected cattle
across 1,600 herds, for the environment-based transmission model.
Provides realistic population structure for demonstrating SISe model
simulations in a cattle disease epidemiology context.

## Usage

``` r
u0_SISe()
```

## Value

A `data.frame` with 1,600 rows (one per herd) and 2 columns:

- S:

  Number of susceptible cattle in the herd

- I:

  Number of infected cattle in the herd (all zero at start)

## Details

This dataset represents initial disease states in a population of 1,600
cattle herds (nodes). Each row represents a single herd (node). The SISe
model extends the SIS model with an environmental compartment for
pathogen shedding, suitable for diseases transmitted through
environmental contamination.

The data contains:

- S:

  Total susceptible cattle in the herd

- I:

  Total infected cattle (initialized to zero)

The herd size distribution reflects realistic heterogeneity observed in
cattle populations, making it suitable for testing environmentally-
mediated transmission dynamics where pathogen survival in the
environment is important.

## See also

[`SISe`](http://stewid.github.io/SimInf/reference/SISe.md) for creating
SISe models with this initial state and
[`events_SISe`](http://stewid.github.io/SimInf/reference/events_SISe.md)
for associated cattle movement and demographic events

## Examples

``` r
## For reproducibility, call the set.seed() function and specify the
## number of threads to use. To use all available threads, remove the
## set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SISe' model with 1600 cattle herds (nodes) and
## initialize it to run over 4*365 days. Add ten infected animals to
## the first herd. Define 'tspan' to record the state of the system at
## weekly time-points. Load scheduled events for the population of
## nodes with births, deaths and between-node movements of
## individuals.
u0 <- u0_SISe()
u0$I[1] <- 10
model <- SISe(u0 = u0,
              tspan = seq(from = 1, to = 4*365, by = 7),
              events = events_SISe(),
              phi = 0,
              upsilon = 1.8e-2,
              gamma = 0.1,
              alpha = 1,
              beta_t1 = 1.0e-1,
              beta_t2 = 1.0e-1,
              beta_t3 = 1.25e-1,
              beta_t4 = 1.25e-1,
              end_t1 = 91,
              end_t2 = 182,
              end_t3 = 273,
              end_t4 = 365,
              epsilon = 0)

## Display the number of cattle affected by each event type per day.
plot(events(model))


## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Plot the median and interquartile range of the number of
## susceptible and infected individuals.
plot(result)


## Plot the trajectory for the first herd.
plot(result, index = 1)


## Summarize the trajectory. The summary includes the number of events
## by event type.
summary(result)
#> Model: SISe
#> Number of nodes: 1600
#> 
#> Transitions
#> -----------
#>  S -> upsilon*phi*S -> I
#>  I -> gamma*I -> S
#> 
#> Global data
#> -----------
#>  Parameter Value
#>  upsilon   0.018
#>  gamma     0.100
#>  alpha     1.000
#>  beta_t1   0.100
#>  beta_t2   0.100
#>  beta_t3   0.125
#>  beta_t4   0.125
#>  epsilon   0.000
#> 
#> Local data
#> ----------
#>  Parameter Value
#>  end_t1     91  
#>  end_t2    182  
#>  end_t3    273  
#>  end_t4    365  
#> 
#> Scheduled events
#> ----------------
#>  Exit: 182535
#>  Enter: 182685
#>  Internal transfer: 0
#>  External transfer: 101472
#> 
#> Network summary
#> ---------------
#>             Min. 1st Qu. Median Mean 3rd Qu. Max.
#>  Indegree:  40.0    57.0   62.0 62.1    68.0 90.0
#>  Outdegree: 36.0    57.0   62.0 62.1    67.0 89.0
#> 
#> Continuous state variables
#> --------------------------
#>       Min. 1st Qu. Median  Mean 3rd Qu.  Max.
#>  phi 0.000   0.000  0.000 0.108   0.000 5.548
#> 
#> Compartments
#> ------------
#>      Min. 1st Qu. Median   Mean 3rd Qu.   Max.
#>  S  18.00  100.00 120.00 122.97  145.00 237.00
#>  I   0.00    0.00   0.00   1.57    0.00 100.00
```
