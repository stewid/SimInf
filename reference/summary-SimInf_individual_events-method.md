# Detailed summary of a `SimInf_individual_events` object

Display a summary of a
[`SimInf_individual_events`](http://stewid.github.io/SimInf/reference/SimInf_individual_events-class.md)
object, including the number of unique individuals, the total number of
event records, and the number of events by event type (exit, enter,
internal transfer, and external transfer).

## Usage

``` r
# S4 method for class 'SimInf_individual_events'
summary(object, ...)
```

## Arguments

- object:

  The
  [`SimInf_individual_events`](http://stewid.github.io/SimInf/reference/SimInf_individual_events-class.md)
  object to summarize.

- ...:

  Additional arguments affecting the summary produced. Currently
  ignored.

## Value

`NULL`, returned invisibly.

## See also

[`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md)
for creating `SimInf_individual_events` objects,
[`SimInf_individual_events`](http://stewid.github.io/SimInf/reference/SimInf_individual_events-class.md)
for the class definition,
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for the node-level event class,
[`show`](http://stewid.github.io/SimInf/reference/show-SimInf_individual_events-method.md)
for a brief summary.
