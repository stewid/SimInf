# Extract the scheduled events from a `SimInf_model` object

Retrieve the `SimInf_events` object containing the schedule of discrete
events (e.g., births, deaths, movements) associated with a
`SimInf_model`. This object holds the timing, location, and type of each
event, as well as the matrices defining how events affect the model
state.

## Usage

``` r
events(object, ...)

# S4 method for class 'SimInf_model'
events(object, ...)
```

## Arguments

- object:

  A `SimInf_model` object.

- ...:

  Additional arguments (currently ignored).

## Value

A
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
object containing the event schedule and associated matrices (`E` and
`N`).

## See also

[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for details on the structure of the returned event object (slots `E`
(select matrix), `N` (shift matrix), `event`, etc.).
[`events_SIR`](http://stewid.github.io/SimInf/reference/events_SIR.md),
[`events_SEIR`](http://stewid.github.io/SimInf/reference/events_SEIR.md),
[`events_SISe3`](http://stewid.github.io/SimInf/reference/events_SISe3.md)
for examples of pre-defined event datasets.
[`mparse`](http://stewid.github.io/SimInf/reference/mparse.md) for
defining custom models with event schedules.
[`run`](http://stewid.github.io/SimInf/reference/run.md) for executing
the simulation with the scheduled events. Vignette `"Scheduled events"`
for a comprehensive tutorial on defining event data, using the
**select** (`E`) and **shift** (`N`) matrices, and simulating complex
demographic and movement processes.

## Examples

``` r
## Create an SIR model with scheduled events.
model <- SIR(
  u0     = u0_SIR(),
  tspan  = 1:(4 * 365),
  events = events_SIR(),
  beta   = 0.16,
  gamma  = 0.077
)

## Extract the events and display a summary.
ev <- events(model)
summary(ev)
#> Number of scheduled events: 466692
#>  - Exit: 182535 (n: min = 1 max = 1 avg = 1.0)
#>  - Enter: 182685 (n: min = 1 max = 1 avg = 1.0)
#>  - Internal transfer: 0
#>  - External transfer: 101472 (n: min = 1 max = 1 avg = 1.0)

## Plot the event schedule over time.
plot(ev)
```
