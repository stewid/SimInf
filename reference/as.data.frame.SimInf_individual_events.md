# Coerce a `SimInf_individual_events` object to a `data.frame`

Convert the cleaned individual-level events stored in a
`SimInf_individual_events` object into a `data.frame`. This function
extracts the individual identifier, event type, time, source node, and
destination node. The resulting `data.frame` has one row per event.

## Usage

``` r
# S3 method for class 'SimInf_individual_events'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `SimInf_individual_events` object.

- ...:

  Additional arguments (currently ignored).

## Value

A `data.frame` with columns:

- `id`: Identifier of the individual (integer or character, depending on
  input).

- `event`: Event type (integer or character, depending on input).

- `time`: Time of the event (numeric or `Date`, depending on input).

- `node`: Source node identifier (integer or character, depending on
  input).

- `dest`: Destination node identifier (integer or character, depending
  on input) (may be `NA`).

## See also

[`SimInf_individual_events`](http://stewid.github.io/SimInf/reference/SimInf_individual_events-class.md)
for the class definition and
[`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md)
for cleaning livestock event data and prepare it for usage in SimInf.
