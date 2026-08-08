# Brief summary of a `SimInf_events` object

Display the number of scheduled events in a
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
object. The count reflects the total number of event entries in the
event schedule, including both internal (E1) and external (E2) events.

## Usage

``` r
# S4 method for class 'SimInf_events'
show(object)
```

## Arguments

- object:

  The
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
  object to display.

## Value

The `object`, returned invisibly.

## See also

[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events.md)
for creating event objects,
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for the class definition,
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for how events are attached to a model.

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

## Brief summary of the events
events(model)
#> Number of scheduled events: 466692
```
