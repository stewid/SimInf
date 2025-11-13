# Transform individual events to node events for a model

In many countries, individual-based livestock data are collected to
enable contact tracing during disease outbreaks. However, the livestock
databases are not always structured in such a way that relevant
information for disease spread simulations is easily retrieved. The aim
of this function is to facilitate cleaning livestock event data and
prepare it for usage in SimInf.

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

  All events that occur after ‘time’ are included. Default is `NULL`
  which means to extract the events after the minimum time-point in the
  `SimInf_individual_events` object.

- target:

  The SimInf model ('SEIR', 'SIR', 'SIS', 'SISe3', 'SISe3_sp', 'SISe',
  or 'SISe_sp') to target the events and u0 for. The default, `NULL`,
  creates events but they might have to be post-processed to fit the
  specific use case.

- age:

  Integer vector with break points in days for the ageing events.

## Value

a `data.frame` with the columns `event`, `time`, `node`, `dest`, `n`,
`proportion`, `select`, and `shift`.

## Details

The individual-based events will be aggregated on node-level. The
`select` value is determined by the event type and age category. If
there is only one age category, i.e., `age=NULL`, then `select=1` for
the enter events, and `select=2` for all other events. If there are two
age categories, then `select=1` for the enter events in the first age
category, and `select=2` for the enter events in the second age
category. Similarly, `select=3` for all other events in the first age
category, and `select=4` for all other events in the first second
category. With three age categories, it works similarly with
`select=1,2,3` for the enter events in each age category, respectively.
And `select=4,5,6` for all other events.

## See also

[`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md).
