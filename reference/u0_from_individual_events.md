# Derive the initial compartment state from individual events

Compute the initial number of individuals in each compartment for each
node based on a set of individual events. The function sums the net
effect of all events occurring before a specified `time` point to
determine the starting state of the simulation.

## Usage

``` r
u0_from_individual_events(events, time = NULL, target = NULL, age = NULL)

# S4 method for class 'SimInf_individual_events'
u0_from_individual_events(events, time = NULL, target = NULL, age = NULL)

# S4 method for class 'data.frame'
u0_from_individual_events(events, time = NULL, target = NULL, age = NULL)
```

## Arguments

- events:

  A `SimInf_individual_events` object containing cleaned individual
  events (e.g., births, deaths, movements) processed by
  [`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md).

- time:

  A numeric scalar specifying the time point at which to calculate the
  initial state. If `NULL` (default), the earliest time point among the
  events is used.

- target:

  A character string specifying the target model type (e.g., `"SIR"`,
  `"SEIR"`, `"SISe3"`). If provided, the function ensures the output
  `u0` includes all required compartments for that model (e.g., adding
  zero-initialized `I`, `R`, or age-specific compartments like `S_1`,
  `S_2`). If `NULL` (default), the output contains only the compartments
  derived from the events (typically `S` or age-stratified `S_*`), which
  may require manual renaming or expansion for specific models.

- age:

  An integer vector of break points (in days) defining age categories.
  The intervals are defined as:

  - Category 1: Age \< `age[1]`

  - Category 2: `age[1]` \<= Age \< `age[2]`

  - ...

  - Last Category: Age \>= `tail(age, 1)`

  If `NULL` (default), all individuals are assigned to a single
  non-age-stratified susceptible compartment (`S`). If provided, the
  output will include columns `S_1`, `S_2`, etc., corresponding to the
  defined age intervals.

## Value

A `data.frame` with one row per node and columns representing the
initial state. The columns include:

- `key`: The original node identifier from the events.

- `node`: A sequential integer node index (1, 2, ...).

- `S_*`: Columns for susceptible individuals (either `S` or
  age-stratified `S_1`, `S_2`, etc.).

- Additional compartments (e.g., `I`, `R`, `E`) if `target` is
  specified, initialized to zero.

## Details

This function accepts two types of input for the `events` argument:

- A `SimInf_individual_events` object (already cleaned by
  [`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md)).

- A raw `data.frame` of events. If a data frame is provided, it is
  automatically cleaned and processed using
  [`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md)
  before the initial state is calculated.

This is particularly useful for initializing models from historical
movement or demographic data, ensuring the simulation starts with the
correct population structure derived from the event log.

## See also

[`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md)
for processing raw event data,
[`u0`](http://stewid.github.io/SimInf/reference/u0.md) for retrieving
the initial state of a model, and `u0<-` for updating the initial state.
