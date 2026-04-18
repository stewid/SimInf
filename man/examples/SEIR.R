## For reproducibility, call the set.seed() function and specify the
## number of threads to use. To use all available threads, remove the
## set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create a 'SEIR' model with 1600 cattle herds (nodes) and initialize
## it to run over 4*365 days. Add ten exposed animals to the first
## herd. Define 'tspan' to record the state of the system at weekly
## time-points. Load scheduled events for the population of nodes with
## births, deaths and between-node movements of individuals.
u0 <- u0_SEIR()
u0$E[1] <- 10
model <- SEIR(
    u0      = u0,
    tspan   = seq(from = 1, to = 4*365, by = 7),
    events  = events_SEIR(),
    beta    = 0.16,
    epsilon = 0.25,
    gamma   = 0.01
)

## Display the number of cattle affected by each event type per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Plot the median and interquartile range of the number of
## susceptible, exposed, infected and recovered individuals.
plot(result)

## Plot the trajectory for the first herd.
plot(result, index = 1)

## Summarize the trajectory. The summary includes the number of events
## by event type.
summary(result)
