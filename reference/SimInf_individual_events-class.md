# Class `SimInf_individual_events`

Storage class for cleaned individual-level event data, such as births,
deaths, and movements. This class serves as an intermediate step in the
data preparation pipeline, holding raw records that will later be
aggregated into the scheduled events
([`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md))
required by a
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
at predefined time-points.

## Details

The typical workflow is:

1.  Collect raw individual event data (e.g., from a database).

2.  Clean and validate it using
    [`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md),
    which returns a `SimInf_individual_events` object.

3.  Aggregate the individual events into node-level scheduled events for
    the model, for example using
    [`u0_from_individual_events`](http://stewid.github.io/SimInf/reference/u0_from_individual_events.md)
    to derive the initial state.

## Slots

- `id`:

  An integer or character vector serving as a unique identifier for each
  individual.

- `event`:

  The type of event. Four types are supported: *exit*, *enter*,
  *internal transfer*, and *external transfer*. These can be specified
  as either a numeric code or a character string:

  - `0` or `"exit"`: Individual leaves the system.

  - `1` or `"enter"`: Individual enters the system.

  - `2` or `"intTrans"`: Individual moves within the same node.

  - `3` or `"extTrans"`: Individual moves to a different node.

- `time`:

  A numeric, character, or `Date` vector indicating when the event
  occurred. Character strings must be coercible to `Date` (e.g.,
  "2023-01-15").

- `node`:

  An integer or character vector identifying the **source** node for the
  event.

- `dest`:

  An integer or character vector identifying the **destination** node.
  For *exit* and *enter* events, this value is typically `NA` or unused.

## See also

[`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md)
for cleaning and processing raw event data into this class format,
[`u0_from_individual_events`](http://stewid.github.io/SimInf/reference/u0_from_individual_events.md)
for deriving the initial state from these events, and
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for the node-level aggregated event format used by the solver.
