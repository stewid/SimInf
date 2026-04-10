# Example event data for the SIS model with cattle herds

Dataset containing 466,692 scheduled events for a population of 1,600
cattle herds over 1,460 days (4 years). Demonstrates how demographic and
movement events affect SIS dynamics in a cattle disease context.

## Usage

``` r
events_SIS()
```

## Value

A `data.frame` with columns:

- event:

  Event type: "exit", "enter", or "extTrans".

- time:

  Day when event occurs (1-1460).

- node:

  Affected herd identifier (1-1600).

- dest:

  Destination herd for external transfer events.

- n:

  Number of cattle affected.

- proportion:

  0\. Not used in this example.

- select:

  Model compartment to affect (see
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)).

- shift:

  0\. Not used in this example.

## Details

The event data contains three types of scheduled events that affect
cattle herds (nodes):

- Exit:

  Deaths or removal of cattle from a herd (n = 182,535). These events
  decrease the population in both susceptible and infected compartments.

- Enter:

  Births or introduction of cattle to a herd (n = 182,685). These events
  add susceptible cattle to herds.

- External transfer:

  Movement of cattle between herds (n = 101,472). These events transfer
  cattle from one herd to another, potentially spreading disease across
  the herd network. Either susceptible or infected animals may be
  transferred.

The `select` column in the returned data frame is mapped to the columns
of the internal select matrix:

- `select = 1` corresponds to **Enter** events, targeting the
  Susceptible (S) compartment.

- `select = 2` corresponds to **Exit** and **External Transfer** events,
  targeting all compartments (S and I).

Events are distributed across all 1,600 herds over the 4-year period.
These are synthetic data generated to illustrate how to incorporate
scheduled events (such as births, deaths, and movements) into a
compartment model in the SimInf framework.

## See also

[`u0_SIS`](http://stewid.github.io/SimInf/reference/u0_SIS.md) for the
corresponding initial cattle population,
[`SIS`](http://stewid.github.io/SimInf/reference/SIS.md) for creating
SIS models with these events, and
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for event structure details

## Examples

``` r
## For reproducibility, call the set.seed() function and specify the
## number of threads to use. To use all available threads, remove the
## set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIS' model with 1600 cattle herds (nodes) and initialize
## it to run over 4*365 days. Add one infected animal to the first
## herd to seed the outbreak. Define 'tspan' to record the state of
## the system at daily time-points. Load scheduled events for the
## population of nodes with births, deaths and between-node movements
## of individuals.
u0 <- u0_SIS()
u0$I[1] <- 1
model <- SIS(u0     = u0,
             tspan  = seq(from = 1, to = 4*365, by = 1),
             events = events_SIS(),
             beta   = 0.16,
             gamma  = 0.01)

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
#> Model: SIS
#> Number of nodes: 1600
#> 
#> Transitions
#> -----------
#>  S -> beta*S*I/(S+I) -> I
#>  I -> gamma*I -> S
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
#>  S   0.0     7.0   10.0  44.2    96.0 218.0
#>  I   0.0     0.0   96.0  80.3   125.0 228.0
```
