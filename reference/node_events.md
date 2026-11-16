# Transform individual events to node events for a model

In many countries, individual-based livestock data are collected to
enable contact tracing during disease outbreaks. However, the
individual-based records are not always directly usable for disease
spread simulations. This function aggregates the cleaned individual
events from
[`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md)
on node-level to prepare them for usage in SimInf models.

## Usage

``` r
node_events(x, time = NULL, target = NULL, age = NULL)

# S4 method for class 'SimInf_individual_events'
node_events(x, time = NULL, target = NULL, age = NULL)
```

## Arguments

- x:

  an individual events object of class `SimInf_individual_events`.

- time:

  All events that occur strictly after ‘time’ are included. Default is
  `NULL` which means to extract the events after the minimum time-point
  in the `SimInf_individual_events` object.

- target:

  The SimInf model ('SEIR', 'SIR', 'SIS', 'SISe3', 'SISe3_sp', 'SISe',
  or 'SISe_sp') to target the events and u0 for. The default, `NULL`,
  creates events but they might have to be post-processed to fit the
  specific use case.

- age:

  Integer vector with break points in days for the ageing events, see
  Details.

## Value

a `data.frame` with the columns `event`, `time`, `node`, `dest`, `n`,
`proportion`, `select`, and `shift`.

## Details

The nodes in the returned events are numbered 1 to N in sorted order of
the unique node identifiers in the individual events. If the individual
events were specified with character event names or dates, the `event`
and `time` columns of the returned events are converted back to
character.

All individuals must have an 'enter' event, otherwise an error is
raised.

The individual-based events are aggregated on node-level. The `select`
value is determined by the event type and age category. If there is only
one age category, i.e., `age=NULL`, then `select=1` for the enter
events, and `select=2` for all other events when the events are targeted
for a model. If there are two age categories, then `select=1` for the
enter events in the first age category, and `select=2` for the enter
events in the second age category. Similarly, `select=3` for all other
events in the first age category, and `select=4` for all other events in
the second age category. With three age categories, it works similarly
with `select=1,2,3` for the enter events in each age category,
respectively. And `select=4,5,6` for all other events.

For the built-in models, the `select` value for the non-enter events is
determined by the design of the target model, such that the events act
on the appropriate compartment.

If `age` contains more than one break point, internal transfer events
for the ageing of individuals between age categories are generated and
included in the returned events. The age categories are determined by
the break points in `age`: individuals with an age below `age[1]` belong
to the first category, individuals with an age from `age[i-1]` up to
`age[i]` belong to category `i`, and individuals with an age of at least
the last break point belong to the last category.

## See also

[`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md),
[`u0_from_individual_events`](http://stewid.github.io/SimInf/reference/u0_from_individual_events.md),
[`events`](http://stewid.github.io/SimInf/reference/events.md).
