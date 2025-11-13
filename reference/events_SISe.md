# Example data to initialize events for the ‘SISe’ model

Example data to initialize scheduled events for a population of 1600
nodes and demonstrate the
[`SISe`](http://stewid.github.io/SimInf/reference/SISe-class.md) model.

## Usage

``` r
events_SISe()
```

## Value

A `data.frame`

## Details

Example data to initialize scheduled events (see
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md))
for a population of 1600 nodes and demonstrate the
[`SISe`](http://stewid.github.io/SimInf/reference/SISe-class.md) model.
The dataset contains 466692 events for 1600 nodes distributed over 4 \*
365 days. The events are divided into three types: ‘Exit’ events remove
individuals from the population (n = 182535), ‘Enter’ events add
individuals to the population (n = 182685), and ‘External transfer’
events move individuals between nodes in the population (n = 101472).
The vignette contains a detailed description of how scheduled events
operate on a model.

## Examples

``` r
## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SISe' model with 1600 nodes and initialize
## it to run over 4*365 days. Add one infected individual
## to the first node.
u0 <- u0_SISe()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 4*365, by = 1)
model <- SISe(u0 = u0, tspan = tspan, events = events_SISe(),
              phi = 0, upsilon = 1.8e-2, gamma = 0.1, alpha = 1,
              beta_t1 = 1.0e-1, beta_t2 = 1.0e-1, beta_t3 = 1.25e-1,
              beta_t4 = 1.25e-1, end_t1 = 91, end_t2 = 182,
              end_t3 = 273, end_t4 = 365, epsilon = 0)

## Display the number of individuals affected by each event type
## per day.
plot(events(model))


## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Summarize the trajectory. The summary includes the number of
## events by event type.
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
#>          Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#>  phi 0.00e+00 0.00e+00 0.00e+00 8.92e-08 0.00e+00 1.88e-02
#> 
#> Compartments
#> ------------
#>        Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#>  S 1.80e+01 1.01e+02 1.22e+02 1.25e+02 1.46e+02 2.37e+02
#>  I 0.00e+00 0.00e+00 0.00e+00 1.28e-06 0.00e+00 1.00e+00
```
