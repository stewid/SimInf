\dontrun{
## For reproducibility, call the set.seed() function and specify the
## number of threads to use. To use all available threads, remove the
## set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIS' model with 1600 cattle herds (nodes) and initialize
## it to run over 4*365 days. Add one infected animal to the first
## herd to seed the outbreak. Define 'tspan' to record the state of
## the system at daily time-points. Load scheduled events for the
## population of nodes with births, deaths and between-node movements
## of individuals.
u0 <- u0_SIS()
u0$I[1] <- 1
model <- SIS(u0     = u0,
             tspan  = seq(from = 1, to = 4*365, by = 1),
             events = events_SIS(),
             beta   = 0.16,
             gamma  = 0.01)

## Display the number of cattle affected by each event type per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Plot the median and interquartile range of the number of
## susceptible and infected individuals.
plot(result)

## Plot the trajectory for the first herd.
plot(result, index = 1)

## Summarize the trajectory. The summary includes the number of events
## by event type.
summary(result)
}
