# Individual events

In many countries, individual-based livestock data are collected to
enable contact tracing during disease outbreaks. However, the livestock
databases are not always structured in such a way that relevant
information for disease spread simulations is easily retrieved. The aim
of this function is to facilitate cleaning livestock event data and
prepare it for usage in SimInf.

## Usage

``` r
individual_events(events)
```

## Arguments

- events:

  a `data.frame` with the columns `id`, `event`, `time`, `node`, and
  `dest` to define the events, see `details`.

## Value

[SimInf_individual_events](http://stewid.github.io/SimInf/reference/SimInf_individual_events-class.md)

## Details

The argument `events` in `individual_events` must be a `data.frame` with
the following columns:

- **id:** an integer or character identifier of the individual.

- **event:** four event types are supported: *exit*, *enter*, *internal
  transfer*, and *external transfer*. When assigning the events, they
  can either be coded as a numerical value or a character string:
  *exit;* `0` or `'exit'`, *enter;* `1` or `'enter'`, *internal
  transfer;* `2` or `'intTrans'`, and *external transfer;* `3` or
  `'extTrans'`.

- **time:** an integer, character, or date (of class `Date`) for when
  the event occured. If it's a character it must be able to coerce to
  `Date`.

- **node:** an integer or character identifier of the source node.

- **dest:** an integer or character identifier of the destination node
  for movement events, and 'dest' will be replaced with `NA` for
  non-movement events.

## See also

[`node_events`](http://stewid.github.io/SimInf/reference/node_events.md).
