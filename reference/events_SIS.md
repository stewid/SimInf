# Example event data for the SIS model with cattle herds

Dataset containing 466,692 scheduled events for a population of 1,600
cattle herds over 1,460 days (4 years). Demonstrates how demographic and
movement events affect SIS dynamics in a cattle disease context.

## Usage

``` r
events_SIS()
```

## Value

A `data.frame` with columns:

- event:

  Event type: "exit", "enter", or "extTrans"

- time:

  Day when event occurs (1-1460)

- node:

  Affected herd identifier (1-1600)

- dest:

  Destination herd for external transfer events

- n:

  Number of cattle affected

- select:

  Model compartment to affect (see
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md))

## Details

The event data contains three types of scheduled events that affect
cattle herds (nodes):

- Exit:

  Deaths or removal of cattle from a herd (n = 182,535). These events
  decrease the population in both susceptible and infected compartments.

- Enter:

  Births or introduction of cattle to a herd (n = 182,685). These events
  add susceptible cattle to herds.

- External transfer:

  Movement of cattle between herds (n = 101,472). These events transfer
  cattle from one herd to another, potentially spreading disease across
  the herd network. Either susceptible or infected animals may be
  transferred.

Events are distributed across all 1,600 herds over the 4-year period,
reflecting realistic patterns of cattle demographic change and
herd-to-herd movement. In SIS dynamics, these events can introduce
disease to previously unaffected herds or remove infected cattle from
the system.

## See also

[`u0_SIS`](http://stewid.github.io/SimInf/reference/u0_SIS.md) for the
corresponding initial cattle population,
[`SIS`](http://stewid.github.io/SimInf/reference/SIS.md) for creating
SIS models with these events, and
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for event structure details
