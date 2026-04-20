# Coerce a `SimInf_events` object to a `data.frame`

Convert the scheduled events stored in a `SimInf_events` object into a
standard `data.frame`. This function extracts the event type, time,
source node, destination node, number of individuals, proportion, and
the specific columns from the select (`E`) and shift (`N`) matrices that
define how each event modifies the compartment state. The resulting
`data.frame` has one row per scheduled event.

## Usage

``` r
# S3 method for class 'SimInf_events'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `SimInf_events` object.

- ...:

  Additional arguments (currently ignored).

## Value

A `data.frame` with columns:

- `event`: Event type (numeric or character, depending on input).

- `time`: Time of the event (numeric or `Date`, depending on input).

- `node`: Source node identifier.

- `dest`: Destination node identifier (may be `NA`).

- `n`: Number of individuals affected.

- `proportion`: Proportion of the population affected (if applicable).

- `select`: The column vector from the select matrix (`E`) that defines
  how the event modifies the compartment state.

- `shift`: The column vector from the shift matrix (`N`) that defines
  how the event modifies the compartment state.

## See also

[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
for the class definition and
[`events`](http://stewid.github.io/SimInf/reference/events.md) for
extracting events from a model.

## Examples

``` r
## Create an 'SIR' model with 1600 cattle herds (nodes) and
## initialize it to run over 4*365 days. Define 'tspan' to record
## the state of the system at daily time-points. Load scheduled
## events for the population of nodes with births, deaths and
## between-node movements of individuals.
model <- SIR(
  u0     = u0_SIR(),
  tspan  = seq(from = 1, to = 4*365, by = 1),
  events = events_SIR(),
  beta   = 0.16,
  gamma  = 0.01
)

## Extract the events from the model and convert to a data frame.
head(as.data.frame(events(model)))
#>   event time node dest n proportion select shift
#> 1  exit    1    6    0 1          0      4     0
#> 2  exit    1   12    0 1          0      4     0
#> 3  exit    1   30    0 1          0      4     0
#> 4  exit    1   46    0 1          0      4     0
#> 5  exit    1   55    0 1          0      4     0
#> 6  exit    1   57    0 1          0      4     0
```
