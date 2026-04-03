# Example event data for the SISe model with cattle herds

Dataset containing 466,692 scheduled events for a population of 1,600
cattle herds over 1,460 days (4 years). Demonstrates how demographic and
movement events affect SISe dynamics in a cattle disease context.

## Usage

``` r
events_SISe()
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
cattle herds:

- Exit:

  Deaths or removal of cattle from a herd (n = 182,535). These events
  decrease the population in susceptible and infected compartments.

- Enter:

  Births or introduction of cattle to a herd (n = 182,685). These events
  add susceptible cattle to herds.

- External transfer:

  Movement of cattle between herds (n = 101,472). These events transfer
  cattle from one herd to another, potentially introducing infected
  animals.

Events are distributed across all 1,600 herds over the 4-year period,
reflecting realistic patterns of cattle demographic change and
herd-to-herd movement.

## See also

[`u0_SISe`](http://stewid.github.io/SimInf/reference/u0_SISe.md) for the
corresponding initial cattle population,
[`SISe`](http://stewid.github.io/SimInf/reference/SISe.md) for creating
SISe models with these events and
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
} # }
```
