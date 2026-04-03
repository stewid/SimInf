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
if (FALSE) { # \dontrun{
## For reproducibility, call the set.seed() function and specify the
## number of threads to use. To use all available threads, remove the
## set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 1600 cattle herds (nodes) and initialize
## it to run over 4*365 days. Add one infected animal to the first
## herd to seed the outbreak. Define 'tspan' to record the state of
## the system at daily time-points. Load scheduled events for the
## population of nodes with births, deaths and between-node movements
## of individuals.
u0 <- u0_SIR()
u0$I[1] <- 1
model <- SIR(u0     = u0,
             tspan  = seq(from = 1, to = 4*365, by = 1),
             events = events_SIR(),
             beta   = 0.16,
             gamma  = 0.01)

## Display the number of cattle affected by each event type per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Plot the median and interquartile range of the number of
## susceptible, infected and recovered individuals.
plot(result)

## Plot the trajectory for the first herd.
plot(result, index = 1)

## Summarize the trajectory. The summary includes the number of events
## by event type.
summary(result)
} # }
```
