# Example event data for the SEIR model with cattle herds

Dataset containing 466,692 scheduled events for a population of 1,600
cattle herds over 1,460 days (4 years). Demonstrates how demographic and
movement events affect SEIR dynamics in a cattle disease context.

## Usage

``` r
events_SEIR()
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
  decrease the population and affect all disease compartments
  proportionally.

- Enter:

  Births or introduction of cattle to a herd (n = 182,685). These events
  add susceptible cattle to herds, increasing overall herd size.

- External transfer:

  Movement of cattle between herds (n = 101,472). These events transfer
  cattle from one herd to another, potentially facilitating between-herd
  disease transmission.

Events are distributed across all 1,600 herds over the 4-year period,
reflecting realistic patterns of cattle demographic change and
herd-to-herd movement. The timing and frequency of events can
significantly influence disease dynamics simulated by the model.

## See also

[`u0_SEIR`](http://stewid.github.io/SimInf/reference/u0_SEIR.md) for the
corresponding initial cattle population,
[`SEIR`](http://stewid.github.io/SimInf/reference/SEIR.md) for creating
SEIR models with these events, and
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

## Create a 'SEIR' model with 1600 cattle herds (nodes) and initialize
## it to run over 4*365 days. Add ten exposed animals to the first
## herd. Define 'tspan' to record the state of the system at weekly
## time-points. Load scheduled events for the population of nodes with
## births, deaths and between-node movements of individuals.
u0 <- u0_SEIR()
u0$E[1] <- 10
model <- SEIR(u0      = u0,
              tspan   = seq(from = 1, to = 4*365, by = 7),
              events  = events_SEIR(),
              beta    = 0.16,
              epsilon = 0.25,
              gamma   = 0.01)

## Display the number of cattle affected by each event type per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Plot the median and interquartile range of the number of
## susceptible, exposed, infected and recovered individuals.
plot(result)

## Plot the trajectory for the first herd.
plot(result, index = 1)

## Summarize the trajectory. The summary includes the number of events
## by event type.
summary(result)
} # }
```
