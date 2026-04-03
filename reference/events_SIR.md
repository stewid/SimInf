# Example Event Data for the SIR Model with Cattle Herds

Dataset containing 466,692 scheduled events for a population of 1,600
cattle herds over 1,460 days (4 years). Demonstrates how demographic and
movement events affect SIR dynamics in a cattle disease context.

## Usage

``` r
events_SIR()
```

## Value

A `data.frame` with columns:

- event:

  Event type: "exit", "enter", or "extTrans"

- time:

  Day when event occurs (1-1460)

- node:

  Affected herd identifier (1-1600)

- dest:

  Destination herd for external transfer events

- n:

  Number of cattle affected

- select:

  Model compartment to affect (see
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md))

## Details

The event data contains three types of scheduled events that affect
cattle herds (nodes):

- Exit:

  Deaths or removal of cattle from a herd (n = 182,535). These events
  decrease the population and remove cattle from the disease system.

- Enter:

  Births or introduction of cattle to a herd (n = 182,685). These events
  add susceptible cattle to herds, increasing potential targets for
  infection.

- External transfer:

  Movement of cattle between herds (n = 101,472). These events transfer
  cattle from one herd to another, potentially spreading disease across
  the herd network.

Events are distributed across all 1,600 herds over the 4-year period,
reflecting realistic patterns of cattle demographic change and
herd-to-herd movement in a livestock production system.

## See also

[`u0_SIR`](http://stewid.github.io/SimInf/reference/u0_SIR.md) for the
corresponding initial cattle population,
[`SIR`](http://stewid.github.io/SimInf/reference/SIR.md) for creating
SIR models with these events, and
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for event structure details

## Examples

``` r
## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 1600 cattle herds and initialize
## it to run over 4*365 days. Add one infected animal to the
## first herd to seed the outbreak.
u0 <- u0_SIR()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 4*365, by = 1)
model <- SIR(u0     = u0,
             tspan  = tspan,
             events = events_SIR(),
             beta   = 0.16,
             gamma  = 0.01)

## Display the number of cattle affected by each event type per day.
plot(events(model))


## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Summarize the trajectory.
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
