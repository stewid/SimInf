# Individual events

In many countries, individual-based livestock data are collected to
enable contact tracing during disease outbreaks. However, the
individual-based records are not always directly usable for disease
spread simulations. The aim of this function is to facilitate cleaning
livestock event data and prepare it for usage in SimInf.

## Usage

``` r
individual_events(events)
```

## Arguments

- events:

  a `data.frame` with the columns `id`, `event`, `time`, `node`, and
  `dest` to define the events, see Details.

## Value

A
[SimInf_individual_events](http://stewid.github.io/SimInf/reference/SimInf_individual_events-class.md)
object with the cleaned events for all individuals with a plausible life
history.

## Details

The argument `events` in `individual_events` must be a `data.frame` with
the following columns:

- **id:** an integer or character identifier of the individual.

- **event:** three event types are supported: *exit*, *enter*, and
  *external transfer*. When assigning the events, they can either be
  coded as a numerical value or a character string: *exit;* `0` or
  `'exit'`, *enter;* `1` or `'enter'`, and *external transfer;* `3` or
  `'extTrans'`.

- **time:** an integer, character, or date (of class `Date`) for when
  the event occurred. If it's a character it must be able to coerce to
  `Date`.

- **node:** an integer or character identifier of the source node.

- **dest:** an integer or character identifier of the destination node
  for movement events, and 'dest' will be replaced with `NA` for
  non-movement events.

The node and dest identifiers must uniquely identify the nodes, e.g.
'XY010', 'XY023', etc. Character identifiers are internally mapped to
consecutive integer node numbers.

For each individual, duplicated records are removed and the function
then searches for the longest plausible sequence of events: an enter
event, followed by zero or more external transfer events, and finally an
exit event, if the individual has one. In a plausible sequence, the
events occur in chronological order and each external transfer departs
from the node where the individual was last registered. All events that
are not part of the longest plausible sequence are discarded, and if no
plausible sequence can be found for an individual, all events for that
individual are dropped. The method is described in Ewerlöf and others
(2025).

Internal transfer events, such as ageing of an individual, are not part
of the event data that this function accepts. Instead, internal transfer
events can be generated when the cleaned events are applied to a model,
see
[`node_events`](http://stewid.github.io/SimInf/reference/node_events.md)
and its `age` argument.

## References

I. R. Ewerlöf, J. Frössling, M. Tråvén, S. Gunnarsson, L. Stengärde, E.
Hurri, and S. Widgren. Exploring structural changes in the Swedish
cattle population and between-holding movements. *Preventive Veterinary
Medicine*, **243**, 106608, 2025.
[doi:10.1016/j.prevetmed.2025.106608](https://doi.org/10.1016/j.prevetmed.2025.106608)

## See also

[`node_events`](http://stewid.github.io/SimInf/reference/node_events.md).

## Examples

``` r
## Example data with two animals; the first has a plausible life
## history, the second does not since the animal is not registered
## at the node of the exit event.
events <- data.frame(
  id = c(1, 1, 1, 1,
         2, 2, 2),
  event = c("enter", "extTrans", "exit", "extTrans",
            "enter", "extTrans", "exit"),
  time = c("2025-03-01", "2025-06-15", "2025-09-01", "2025-12-01",
           "2025-04-01", "2025-07-01", "2025-10-01"),
  node = c("N1", "N1", "N2", "N2",
           "N3", "N3", "N5"),
  dest = c("-", "N2", "-", "N3",
           "-", "N4", "-")
)

## Keep the longest plausible life history for each animal: the
## movement after the exit event of the first animal is discarded,
## and all events of the second animal are dropped. Note that the
## '-' in 'dest' is replaced with NA for non-movement events.
as.data.frame(individual_events(events))
#>   id    event       time node dest
#> 1  1    enter 2025-03-01   N1 <NA>
#> 2  1 extTrans 2025-06-15   N1   N2
#> 3  1     exit 2025-09-01   N2 <NA>
```
