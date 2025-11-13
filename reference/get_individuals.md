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

  the time-point for the lookup of individuals. Default is `NULL` which
  means to extract the individuals at the minimum time-point in the
  events object `x`.

## Value

a `data.frame` with the columns `id`, `node`, and `age`.
