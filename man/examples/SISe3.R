## For reproducibility, call the set.seed() function and specify the
## number of threads to use. To use all available threads, remove the
## set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SISe3' model with 1600 cattle herds (nodes) stratified
## by age, initialize it to run over 4*365 days and record data at
## weekly time-points. Add ten infected animals to age category 1 in
## the first herd to seed the outbreak.  Define 'tspan' to record the
## state of the system at weekly time-points. Load scheduled events
## events for the population of nodes with births, deaths and
## between-node movements of individuals.
u0 <- u0_SISe3
u0$I_1[1] <- 10
model <- SISe3(
    u0 = u0,
    tspan = seq(from = 1, to = 4*365, by = 7),
    events = events_SISe3,
    phi = rep(0, nrow(u0)),
    upsilon_1 = 1.8e-2,
    upsilon_2 = 1.8e-2,
    upsilon_3 = 1.8e-2,
    gamma_1 = 0.1,
    gamma_2 = 0.1,
    gamma_3 = 0.1,
    alpha = 1,
    beta_t1 = 1.0e-1,
    beta_t2 = 1.0e-1,
    beta_t3 = 1.25e-1,
    beta_t4 = 1.25e-1,
    end_t1 = 91,
    end_t2 = 182,
    end_t3 = 273,
    end_t4 = 365,
    epsilon = 0
)

## Display the number of cattle affected by each event type per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Plot the median and interquartile range of the number of
## susceptible and infected individuals.
plot(result)

## Plot the proportion of nodes with at least one infected individual.
plot(result, I_1 + I_2 + I_3 ~ ., level = 2, type = "l")

## Plot the trajectory for the first herd.
plot(result, index = 1)

## Summarize the trajectory. The summary includes the number of events
## by event type.
summary(result)
