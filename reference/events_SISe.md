# Example data to initialize events for the ‘SISe’ model

Example data to initialize scheduled events for a population of 1600
nodes and demonstrate the
[`SISe`](http://stewid.github.io/SimInf/reference/SISe-class.md) model.

## Usage

``` r
events_SISe()
```

## Value

A `data.frame`

## Details

Example data to initialize scheduled events (see
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md))
for a population of 1600 nodes and demonstrate the
[`SISe`](http://stewid.github.io/SimInf/reference/SISe-class.md) model.
The dataset contains 466692 events for 1600 nodes distributed over 4 \*
365 days. The events are divided into three types: ‘Exit’ events remove
individuals from the population (n = 182535), ‘Enter’ events add
individuals to the population (n = 182685), and ‘External transfer’
events move individuals between nodes in the population (n = 101472).
The vignette contains a detailed description of how scheduled events
operate on a model.

## Examples

``` r
if (FALSE) { # \dontrun{
## For reproducibility, call the set.seed() function and specify the
## number of threads to use. To use all available threads, remove the
## set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SISe' model with 1600 cattle herds (nodes) and
## initialize it to run over 4*365 days. Add ten infected animals to
## the first herd. Define 'tspan' to record the state of the system at
## weekly time-points. Load scheduled events for the population of
## nodes with births, deaths and between-node movements of
## individuals.
u0 <- u0_SISe()
u0$I[1] <- 10
model <- SISe(u0 = u0,
              tspan = seq(from = 1, to = 4*365, by = 7),
              events = events_SISe(),
              phi = 0,
              upsilon = 1.8e-2,
              gamma = 0.1,
              alpha = 1,
              beta_t1 = 1.0e-1,
              beta_t2 = 1.0e-1,
              beta_t3 = 1.25e-1,
              beta_t4 = 1.25e-1,
              end_t1 = 91,
              end_t2 = 182,
              end_t3 = 273,
              end_t4 = 365,
              epsilon = 0)

## Display the number of cattle affected by each event type per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Plot the median and interquartile range of the number of
## susceptible, infected and recovered individuals.
plot(result)

## Plot the trajectory for the first herd.
plot(result, index = 1)

## Summarize the trajectory. The summary includes the number of events
## by event type.
summary(result)
} # }
```
