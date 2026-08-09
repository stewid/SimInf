# Brief summary of a `SimInf_individual_events` object

Display the number of unique individuals and the total number of event
records in a
[`SimInf_individual_events`](http://stewid.github.io/SimInf/reference/SimInf_individual_events-class.md)
object. The individual count reflects the distinct `id` values, while
the event count represents all recorded events (exit, enter, internal
transfer, and external transfer).

## Usage

``` r
# S4 method for class 'SimInf_individual_events'
show(object)
```

## Arguments

- object:

  The
  [`SimInf_individual_events`](http://stewid.github.io/SimInf/reference/SimInf_individual_events-class.md)
  object to display.

## Value

The `object`, returned invisibly.

## See also

[`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md)
for creating `SimInf_individual_events` objects,
[`SimInf_individual_events`](http://stewid.github.io/SimInf/reference/SimInf_individual_events-class.md)
for the class definition,
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for the node-level event class.
