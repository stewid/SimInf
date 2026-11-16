# Extract individuals from `SimInf_individual_events`

Lookup individuals, in which node they are located and their age at a
specified time-point.

## Usage

``` r
get_individuals(x, time = NULL)

# S4 method for class 'SimInf_individual_events'
get_individuals(x, time = NULL)
```

## Arguments

- x:

  an individual events object of class `SimInf_individual_events`.

- time:

  the time-point for the lookup of individuals, either an integer or, if
  the events were specified with dates, a date. Default is `NULL` which
  means to extract the individuals at the minimum time-point in the
  events object `x`.

## Value

a `data.frame` with the columns `id`, `node`, and `age`.

## Details

All individuals must have an 'enter' event, otherwise an error is
raised. An individual is included in the result if its 'enter' event
occurs at or before the time-point, and it has no 'exit' event at or
before the time-point. Individuals that have not yet entered are
omitted, and so are individuals that have exited.

The node is the node where the individual was last registered; if the
last event is an external transfer, the node is the destination of that
movement. Node identifiers are returned as specified in the events, e.g.
'XY010'.

The age is the difference between the time-point and the 'enter' event
of the individual. If the events were specified with dates, the age is
in days, otherwise it is in the same unit as the time-points.

The rows in the returned data frame are ordered by `id`.

## See also

[`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md),
[`node_events`](http://stewid.github.io/SimInf/reference/node_events.md),
[`u0_from_individual_events`](http://stewid.github.io/SimInf/reference/u0_from_individual_events.md).
