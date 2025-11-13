# Example data to initialize events for the ‘SIR’ model

Example data to initialize scheduled events for a population of 1600
nodes and demonstrate the
[`SIR`](http://stewid.github.io/SimInf/reference/SIR-class.md) model.

## Usage

``` r
events_SIR()
```

## Value

A `data.frame`

## Details

Example data to initialize scheduled events (see
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md))
for a population of 1600 nodes and demonstrate the
[`SIR`](http://stewid.github.io/SimInf/reference/SIR-class.md) model.
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

## Create an 'SIR' model with 1600 nodes and initialize
## it to run over 4*365 days. Add one infected individual
## to the first node.
u0 <- u0_SIR()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 4*365, by = 1)
model <- SIR(u0     = u0,
             tspan  = tspan,
             events = events_SIR(),
             beta   = 0.16,
             gamma  = 0.01)

## Display the number of individuals affected by each event type
## per day.
plot(events(model))


## Run the model to generate a single stochastic trajectory.
result <- run(model)
plot(result)


## Summarize the trajectory. The summary includes the number of
## events by event type.
summary(result)
#> Model: SIR
#> Number of nodes: 1600
#> 
#> Transitions
#> -----------
#>  S -> beta*S*I/(S+I+R) -> I
#>  I -> gamma*I -> R
#> 
#> Global data
#> -----------
#>  Number of parameters without a name: 0
#>  - None
#> 
#> Local data
#> ----------
#>  Parameter Value
#>  beta      0.16 
#>  gamma     0.01 
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
#> Compartments
#> ------------
#>     Min. 1st Qu. Median  Mean 3rd Qu.  Max.
#>  S   0.0     5.0   13.0  55.6   112.0 219.0
#>  I   0.0     0.0    4.0  10.9    11.0 168.0
#>  R   0.0     0.0   62.0  58.0   105.0 221.0
```
