# Detailed summary of a `SimInf_events` object

Display a summary of the scheduled events in a
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
object. The output includes the total number of events and summary
statistics (count, minimum, maximum, and average number of individuals
affected by each event).

## Usage

``` r
# S4 method for class 'SimInf_events'
summary(object, ...)
```

## Arguments

- object:

  The
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
  object to summarize.

- ...:

  Additional arguments affecting the summary produced. Currently
  ignored.

## Value

None (invisible `NULL`).

## Details

Four event types are reported:

- **Exit**: Events where individuals are removed from a node (e.g.,
  deaths).

- **Enter**: Events where individuals are added to a node (e.g.,
  births).

- **Internal transfer**: Events moving individuals between compartments
  within a node (e.g., vaccination, ageing).

- **External transfer**: Events moving individuals between different
  nodes (e.g., livestock movements).

## See also

[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events.md)
for creating event objects,
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for the class definition,
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for how events are attached to a model,
[`events`](http://stewid.github.io/SimInf/reference/events.md) for
extracting events from a model.

## Examples

``` r
## Create an 'SIR' model with 1600 cattle herds (nodes) and
## scheduled events for the population of nodes with births,
## deaths and between-node movements of individuals. Define
## 'tspan' to record the state of the system at daily time
## points over 4*365 days.
model <- SIR(
  u0     = u0_SIR(),
  tspan  = 1:(4 * 365),
  events = events_SIR(),
  beta   = 0.16,
  gamma  = 0.077
)

## Detailed summary of the events
summary(events(model))
#> Number of scheduled events: 466692
#>  - Exit: 182535 (n: min = 1 max = 1 avg = 1.0)
#>  - Enter: 182685 (n: min = 1 max = 1 avg = 1.0)
#>  - Internal transfer: 0
#>  - External transfer: 101472 (n: min = 1 max = 1 avg = 1.0)
```
