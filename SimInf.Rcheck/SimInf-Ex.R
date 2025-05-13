pkgname <- "SimInf"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SimInf')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("C_code")
### * C_code

flush(stderr()); flush(stdout())

### Name: C_code
### Title: Extract the C code from a 'SimInf_model' object
### Aliases: C_code

### ** Examples

## Use the model parser to create a 'SimInf_model' object that
## expresses an SIR model, where 'b' is the transmission rate and
## 'g' is the recovery rate.
model <- mparse(transitions = c("S -> b*S*I/(S+I+R) -> I", "I -> g*I -> R"),
                compartments = c("S", "I", "R"),
                gdata = c(b = 0.16, g = 0.077),
                u0 = data.frame(S = 99, I = 1, R = 0),
                tspan = 1:10)

## View the C code.
C_code(model)



cleanEx()
nameEx("SEIR")
### * SEIR

flush(stderr()); flush(stdout())

### Name: SEIR
### Title: Create an SEIR model
### Aliases: SEIR

### ** Examples

## Create a SEIR model object.
model <- SEIR(u0 = data.frame(S = 99, E = 0, I = 1, R = 0),
              tspan = 1:100,
              beta = 0.16,
              epsilon = 0.25,
              gamma = 0.077)

## Run the SEIR model and plot the result.
set.seed(3)
result <- run(model)
plot(result)



cleanEx()
nameEx("SIR-class")
### * SIR-class

flush(stderr()); flush(stdout())

### Name: SIR-class
### Title: Definition of the SIR model
### Aliases: SIR-class

### ** Examples

## Create an SIR model object.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the SIR model and plot the result.
set.seed(22)
result <- run(model)
plot(result)



cleanEx()
nameEx("SIR")
### * SIR

flush(stderr()); flush(stdout())

### Name: SIR
### Title: Create an SIR model
### Aliases: SIR

### ** Examples

## Create an SIR model object.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the SIR model and plot the result.
set.seed(22)
result <- run(model)
plot(result)



cleanEx()
nameEx("SIS-class")
### * SIS-class

flush(stderr()); flush(stdout())

### Name: SIS-class
### Title: Definition of the SIS model
### Aliases: SIS-class

### ** Examples

## Create an SIS model object.
model <- SIS(u0 = data.frame(S = 99, I = 1),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the SIS model and plot the result.
set.seed(22)
result <- run(model)
plot(result)



cleanEx()
nameEx("SIS")
### * SIS

flush(stderr()); flush(stdout())

### Name: SIS
### Title: Create an SIS model
### Aliases: SIS

### ** Examples

## Create an SIS model object.
model <- SIS(u0 = data.frame(S = 99, I = 1),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the SIS model and plot the result.
set.seed(22)
result <- run(model)
plot(result)



cleanEx()
nameEx("SimInf_events")
### * SimInf_events

flush(stderr()); flush(stdout())

### Name: SimInf_events
### Title: Create a 'SimInf_events' object
### Aliases: SimInf_events

### ** Examples

## Let us illustrate how movement events can be used to transfer
## individuals from one node to another.  Use the built-in SIR
## model and start with 2 nodes where all individuals are in the
## first node (100 per compartment).
u0 <- data.frame(S = c(100, 0), I = c(100, 0), R = c(100, 0))

## Then create 300 movement events to transfer all individuals,
## one per day, from the first node to the second node. Use the
## fourth column in the select matrix where all compartments
## can be sampled with equal weight.
events <- data.frame(event      = rep("extTrans", 300),
                     time       = 1:300,
                     node       = 1,
                     dest       = 2,
                     n          = 1,
                     proportion = 0,
                     select     = 4,
                     shift      = 0)

## Create an SIR model without disease transmission to
## demonstrate the events.
model <- SIR(u0      = u0,
             tspan  = 1:300,
             events = events,
             beta   = 0,
             gamma  = 0)

## Run the model and plot the number of individuals in
## the second node.  As can be seen in the figure, all
## indivuduals have been moved to the second node when
## t = 300.
plot(run(model), index = 1:2, range = FALSE)

## Let us now double the weight to sample from the 'I'
## compartment and rerun the model.
model@events@E[2, 4] <- 2
plot(run(model), index = 1:2, range = FALSE)

## And much larger weight to sample from the I compartment.
model@events@E[2, 4] <- 10
plot(run(model), index = 1:2, range = FALSE)

## Increase the weight for the R compartment.
model@events@E[3, 4] <- 4
plot(run(model), index = 1:2, range = FALSE)



cleanEx()
nameEx("abc")
### * abc

flush(stderr()); flush(stdout())

### Name: abc
### Title: Approximate Bayesian computation
### Aliases: abc abc,SimInf_model-method

### ** Examples

## Not run: 
##D ## Let us consider an SIR model in a closed population with N = 100
##D ## individuals of whom one is initially infectious and the rest are
##D ## susceptible. First, generate one realisation (with a specified
##D ## seed) from the model with known parameters \code{beta = 0.16} and
##D ## \code{gamma = 0.077}. Then, use \code{abc} to infer the (known)
##D ## parameters from the simulated data.
##D model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##D              tspan = 1:100,
##D              beta = 0.16,
##D              gamma = 0.077)
##D 
##D ## Run the SIR model and plot the number of infectious.
##D set.seed(22)
##D infectious <- trajectory(run(model), "I")$I
##D plot(infectious, type = "s")
##D 
##D ## The distance function to accept or reject a proposal. Each node
##D ## in the simulated trajectory (contained in the 'result' object)
##D ## represents one proposal.
##D distance <- function(result, ...) {
##D     ## Extract the time-series of infectious in each node as a
##D     ## data.frame.
##D     sim <- trajectory(result, "I")
##D 
##D     ## Split the 'sim' data.frame by node and calculate the sum of the
##D     ## squared distance at each time-point for each node.
##D     dist <- tapply(sim$I, sim$node, function(sim_infectious) {
##D         sum((infectious - sim_infectious)^2)
##D     })
##D 
##D     ## Return the distance for each node. Each proposal will be
##D     ## accepted or rejected depending on if the distance is less than
##D     ## the tolerance for the current generation.
##D     dist
##D }
##D 
##D ## Fit the model parameters using ABC-SMC and adaptive tolerance
##D ## selection. The priors for the parameters are specified using a
##D ## formula notation. Here we use a uniform distribtion for each
##D ## parameter with lower bound = 0 and upper bound = 1. Note that we
##D ## use a low number particles here to keep the run-time of the example
##D ## short. In practice you would want to use many more to ensure better
##D ## approximations.
##D fit <- abc(model = model,
##D            priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
##D            n_particles = 100,
##D            n_init = 1000,
##D            distance = distance,
##D            verbose = TRUE)
##D 
##D ## Print a brief summary.
##D fit
##D 
##D ## Display the ABC posterior distribution.
##D plot(fit)
## End(Not run)



cleanEx()
nameEx("boxplot-SimInf_model-method")
### * boxplot-SimInf_model-method

flush(stderr()); flush(stdout())

### Name: boxplot,SimInf_model-method
### Title: Box plot of number of individuals in each compartment
### Aliases: boxplot,SimInf_model-method

### ** Examples

## Create an 'SIR' model with 10 nodes and initialise
## it with 99 susceptible individuals and one infected
## individual. Let the model run over 100 days.
model <- SIR(u0 = data.frame(S = rep(99, 10),
                             I = rep(1, 10),
                             R = rep(0, 10)),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the model and save the result.
result <- run(model)

## Create a boxplot that includes all compartments in all nodes.
boxplot(result)

## Create a boxplot that includes the S and I compartments in
## nodes 1 and 2.
boxplot(result, ~S+I, 1:2)



cleanEx()
nameEx("distance_matrix")
### * distance_matrix

flush(stderr()); flush(stdout())

### Name: distance_matrix
### Title: Create a distance matrix between nodes for spatial models
### Aliases: distance_matrix

### ** Examples

## Generate a grid 10 x 10 and place one node in each cell
## separated by 100m.
nodes <- expand.grid(x = (0:9) * 100, y = (0:9) * 100)
plot(y ~ x, nodes)

## Define the cutoff to only include neighbors within 300m.
d <- distance_matrix(x = nodes$x, y = nodes$y, cutoff = 300)

## View the first 10 rows and columns in the distance matrix
d[1:10, 1:10]



cleanEx()
nameEx("edge_properties_to_matrix")
### * edge_properties_to_matrix

flush(stderr()); flush(stdout())

### Name: edge_properties_to_matrix
### Title: Convert an edge list with properties to a matrix
### Aliases: edge_properties_to_matrix

### ** Examples

## Let us consider the following edge properties.
edges <- data.frame(
    from  = c(  2,    3,     4,  1,   4,    5,   1,   3,   1,   3),
    to    = c(  1,    1,     1,  2,   3,    3,   4,   4,   5,   5),
    rate  = c(0.2, 0.01,  0.79,  1, 0.2, 0.05, 0.2, 0.8, 0.2, 0.8),
    count = c(  5,    5,     5, 50,  10,   10,   5,   5,   5,   5))

## Converting the edge properties into a matrix
edge_properties_to_matrix(edges, 6)

## Gives the following output. The first column contains first the
## properties for the edge from = 2 --> to = 1, where the first
## row is the zero-based index of from, i.e., 1. The second row
## contains the rate=0.2 and the third row count=5. On the fourth
## row starts the next sequence with the values in the second row
## in the edges data.frame. The stop value in the first column is
## on row 10. As can be seen in column 6, there are no edge
## properties for node=6.
##        [,1] [,2]  [,3] [,4] [,5] [,6]
##  [1,]  1.00    0  3.00  0.0  0.0   -1
##  [2,]  0.20    1  0.20  0.2  0.2  NaN
##  [3,]  5.00   50 10.00  5.0  5.0  NaN
##  [4,]  2.00   -1  4.00  2.0  2.0  NaN
##  [5,]  0.01  NaN  0.05  0.8  0.8  NaN
##  [6,]  5.00  NaN 10.00  5.0  5.0  NaN
##  [7,]  3.00  NaN -1.00 -1.0 -1.0  NaN
##  [8,]  0.79  NaN   NaN  NaN  NaN  NaN
##  [9,]  5.00  NaN   NaN  NaN  NaN  NaN
## [10,] -1.00  NaN   NaN  NaN  NaN  NaN



cleanEx()
nameEx("events")
### * events

flush(stderr()); flush(stdout())

### Name: events
### Title: Extract the events from a 'SimInf_model' object
### Aliases: events events,SimInf_model-method

### ** Examples

## Create an SIR model that includes scheduled events.
model <- SIR(u0     = u0_SIR(),
             tspan  = 1:(4 * 365),
             events = events_SIR(),
             beta   = 0.16,
             gamma  = 0.077)

## Extract the scheduled events from the model and display summary
summary(events(model))

## Extract the scheduled events from the model and plot them
plot(events(model))



cleanEx()
nameEx("events_SEIR")
### * events_SEIR

flush(stderr()); flush(stdout())

### Name: events_SEIR
### Title: Example data to initialize events for the 'SEIR' model
### Aliases: events_SEIR

### ** Examples

## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SEIR' model with 1600 nodes and initialize
## it to run over 4*365 days. Add one infected individual
## to the first node.
u0 <- u0_SEIR()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 4*365, by = 1)
model <- SEIR(u0      = u0,
              tspan   = tspan,
              events  = events_SEIR(),
              beta    = 0.16,
              epsilon = 0.25,
              gamma   = 0.01)

## Display the number of individuals affected by each event type
## per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)
plot(result)

## Summarize the trajectory. The summary includes the number of
## events by event type.
summary(result)



cleanEx()
nameEx("events_SIR")
### * events_SIR

flush(stderr()); flush(stdout())

### Name: events_SIR
### Title: Example data to initialize events for the 'SIR' model
### Aliases: events_SIR

### ** Examples

## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 1600 nodes and initialize
## it to run over 4*365 days. Add one infected individual
## to the first node.
u0 <- u0_SIR()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 4*365, by = 1)
model <- SIR(u0     = u0,
             tspan  = tspan,
             events = events_SIR(),
             beta   = 0.16,
             gamma  = 0.01)

## Display the number of individuals affected by each event type
## per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)
plot(result)

## Summarize the trajectory. The summary includes the number of
## events by event type.
summary(result)



cleanEx()
nameEx("events_SIS")
### * events_SIS

flush(stderr()); flush(stdout())

### Name: events_SIS
### Title: Example data to initialize events for the 'SIS' model
### Aliases: events_SIS

### ** Examples

## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIS' model with 1600 nodes and initialize
## it to run over 4*365 days. Add one infected individual
## to the first node.
u0 <- u0_SIS()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 4*365, by = 1)
model <- SIS(u0     = u0,
             tspan  = tspan,
             events = events_SIS(),
             beta   = 0.16,
             gamma  = 0.01)

## Display the number of individuals affected by each event type
## per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)
plot(result)

## Summarize the trajectory. The summary includes the number of
## events by event type.
summary(result)



cleanEx()
nameEx("events_SISe")
### * events_SISe

flush(stderr()); flush(stdout())

### Name: events_SISe
### Title: Example data to initialize events for the 'SISe' model
### Aliases: events_SISe

### ** Examples

## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SISe' model with 1600 nodes and initialize
## it to run over 4*365 days. Add one infected individual
## to the first node.
u0 <- u0_SISe()
u0$I[1] <- 1
tspan <- seq(from = 1, to = 4*365, by = 1)
model <- SISe(u0 = u0, tspan = tspan, events = events_SISe(),
              phi = 0, upsilon = 1.8e-2, gamma = 0.1, alpha = 1,
              beta_t1 = 1.0e-1, beta_t2 = 1.0e-1, beta_t3 = 1.25e-1,
              beta_t4 = 1.25e-1, end_t1 = 91, end_t2 = 182,
              end_t3 = 273, end_t4 = 365, epsilon = 0)

## Display the number of individuals affected by each event type
## per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Summarize the trajectory. The summary includes the number of
## events by event type.
summary(result)



cleanEx()
nameEx("events_SISe3")
### * events_SISe3

flush(stderr()); flush(stdout())

### Name: events_SISe3
### Title: Example data to initialize events for the 'SISe3' model
### Aliases: events_SISe3
### Keywords: dataset

### ** Examples

## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SISe3' model with 1600 nodes and initialize
## it to run over 4*365 days. Add one infected individual
## to the first node.
data("u0_SISe3", package = "SimInf")
data("events_SISe3", package = "SimInf")
u0_SISe3$I_1[1] <- 1
tspan <- seq(from = 1, to = 4*365, by = 1)
model <- SISe3(u0 = u0_SISe3, tspan = tspan, events = events_SISe3,
               phi = rep(0, nrow(u0_SISe3)), upsilon_1 = 1.8e-2,
               upsilon_2 = 1.8e-2, upsilon_3 = 1.8e-2,
               gamma_1 = 0.1, gamma_2 = 0.1, gamma_3 = 0.1,
               alpha = 1, beta_t1 = 1.0e-1, beta_t2 = 1.0e-1,
               beta_t3 = 1.25e-1, beta_t4 = 1.25e-1, end_t1 = 91,
               end_t2 = 182, end_t3 = 273, end_t4 = 365, epsilon = 0)

## Display the number of individuals affected by each event type
## per day.
plot(events(model))

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Summarize the trajectory. The summary includes the number of
## events by event type.
summary(result)



cleanEx()
nameEx("gdata-set")
### * gdata-set

flush(stderr()); flush(stdout())

### Name: gdata<-
### Title: Set a global data parameter for a 'SimInf_model' object
### Aliases: gdata<- gdata<-,SimInf_model-method

### ** Examples

## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Set 'beta' to a new value
gdata(model, "beta") <- 2

## Extract the global data vector that is common to all nodes
gdata(model)



cleanEx()
nameEx("gdata")
### * gdata

flush(stderr()); flush(stdout())

### Name: gdata
### Title: Extract global data from a 'SimInf_model' object
### Aliases: gdata gdata,SimInf_model-method

### ** Examples

## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Set 'beta' to a new value
gdata(model, "beta") <- 2

## Extract the global data vector that is common to all nodes
gdata(model)



cleanEx()
nameEx("indegree")
### * indegree

flush(stderr()); flush(stdout())

### Name: indegree
### Title: Determine in-degree for each node in a model
### Aliases: indegree

### ** Examples

## Create an 'SIR' model with 1600 nodes and initialize
## it with example data.
model <- SIR(u0 = u0_SIR(), tspan = 1:1460, events = events_SIR(),
             beta   = 0.16, gamma  = 0.077)

## Display indegree for each node in the model.
plot(indegree(model))



cleanEx()
nameEx("ldata")
### * ldata

flush(stderr()); flush(stdout())

### Name: ldata
### Title: Extract local data from a node
### Aliases: ldata ldata,SimInf_model-method

### ** Examples

## Create an 'SISe' model with 1600 nodes.
model <- SISe(u0 = u0_SISe(), tspan = 1:100, events = events_SISe(),
              phi = 0, upsilon = 1.8e-2, gamma = 0.1, alpha = 1,
              beta_t1 = 1.0e-1, beta_t2 = 1.0e-1, beta_t3 = 1.25e-1,
              beta_t4 = 1.25e-1, end_t1 = c(91, 101), end_t2 = c(182, 185),
              end_t3 = c(273, 275), end_t4 = c(365, 360), epsilon = 0)

## Display local data from the first two nodes.
ldata(model, node = 1)
ldata(model, node = 2)



cleanEx()
nameEx("mparse")
### * mparse

flush(stderr()); flush(stdout())

### Name: mparse
### Title: Model parser to define new models to run in 'SimInf'
### Aliases: mparse

### ** Examples

## Not run: 
##D ## Use the model parser to create a 'SimInf_model' object that
##D ## expresses the SIR model, where 'beta' is the transmission rate
##D ## and 'gamma' is the recovery rate.
##D model  <- mparse(transitions = c("S -> beta*S*I/N -> I",
##D                                  "I -> gamma*I -> R",
##D                                  "N <- S+I+R"),
##D                  compartments = c("S", "I", "R"),
##D                  gdata = c(beta = 0.16, gamma = 0.077),
##D                  u0 = data.frame(S = 100, I = 1, R = 0),
##D                  tspan = 1:100)
##D 
##D ## Run and plot the result
##D set.seed(22)
##D result <- run(model)
##D plot(result)
## End(Not run)



cleanEx()
nameEx("n_compartments")
### * n_compartments

flush(stderr()); flush(stdout())

### Name: n_compartments
### Title: Determine the number of compartments in a model
### Aliases: n_compartments n_compartments,SimInf_model-method

### ** Examples

## Create an 'SIR' model with 100 nodes, with 99 susceptible,
## 1 infected and 0 recovered in each node.
u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Display the number of compartments in the model.
n_compartments(model)



cleanEx()
nameEx("n_nodes")
### * n_nodes

flush(stderr()); flush(stdout())

### Name: n_nodes
### Title: Determine the number of nodes in a model
### Aliases: n_nodes n_nodes,SimInf_model-method
###   n_nodes,SimInf_pfilter-method n_nodes,SimInf_pmcmc-method

### ** Examples

## Create an 'SIR' model with 100 nodes, with 99 susceptible,
## 1 infected and 0 recovered in each node.
u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Display the number of nodes in the model.
n_nodes(model)



cleanEx()
nameEx("n_replicates")
### * n_replicates

flush(stderr()); flush(stdout())

### Name: n_replicates
### Title: Determine the number of replicates in a model
### Aliases: n_replicates n_replicates,SimInf_model-method

### ** Examples

## Create an 'SIR' model with 100 nodes, with 99 susceptible,
## 1 infected and 0 recovered in each node.
u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Display the number of replicates in the model.
n_replicates(model)



cleanEx()
nameEx("nodes")
### * nodes

flush(stderr()); flush(stdout())

### Name: nodes
### Title: Example data with spatial distribution of nodes
### Aliases: nodes
### Keywords: dataset

### ** Examples

## Not run: 
##D ## For reproducibility, call the set.seed() function and specify
##D ## the number of threads to use. To use all available threads,
##D ## remove the set_num_threads() call.
##D set.seed(123)
##D set_num_threads(1)
##D 
##D ## Create an 'SIR' model with 1600 nodes and initialize
##D ## it to run over 4*365 days. Add one infected individual
##D ## to the first node.
##D u0 <- u0_SIR()
##D u0$I[1] <- 1
##D tspan <- seq(from = 1, to = 4*365, by = 1)
##D model <- SIR(u0     = u0,
##D              tspan  = tspan,
##D              events = events_SIR(),
##D              beta   = 0.16,
##D              gamma  = 0.077)
##D 
##D ## Run the model to generate a single stochastic trajectory.
##D result <- run(model)
##D 
##D ## Determine nodes with one or more infected individuals in the
##D ## trajectory. Extract the 'I' compartment and check for any
##D ## infected individuals in each node.
##D infected <- colSums(trajectory(result, ~ I, format = "matrix")) > 0
##D 
##D ## Display infected nodes in 'blue' and non-infected nodes in 'yellow'.
##D data("nodes", package = "SimInf")
##D col <- ifelse(infected, "blue", "yellow")
##D plot(y ~ x, nodes, col = col, pch = 20, cex = 2)
## End(Not run)



cleanEx()
nameEx("outdegree")
### * outdegree

flush(stderr()); flush(stdout())

### Name: outdegree
### Title: Determine out-degree for each node in a model
### Aliases: outdegree

### ** Examples

## Create an 'SIR' model with 1600 nodes and initialize
## it with example data.
model <- SIR(u0 = u0_SIR(), tspan = 1:1460, events = events_SIR(),
             beta   = 0.16, gamma  = 0.077)

## Display outdegree for each node in the model.
plot(outdegree(model))



cleanEx()
nameEx("pairs-SimInf_model-method")
### * pairs-SimInf_model-method

flush(stderr()); flush(stdout())

### Name: pairs,SimInf_model-method
### Title: Scatterplot of number of individuals in each compartment
### Aliases: pairs,SimInf_model-method

### ** Examples

## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 10 nodes and initialise
## it with 99 susceptible individuals and one infected
## individual. Let the model run over 100 days.
model <- SIR(u0 = data.frame(S = rep(99, 10),
                             I = rep(1, 10),
                             R = rep(0, 10)),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the model and save the result.
result <- run(model)

## Create a scatter plot that includes all compartments in all
## nodes.
pairs(result)

## Create a scatter plot that includes the S and I compartments in
## nodes 1 and 2.
pairs(result, ~S+I, 1:2)



cleanEx()
nameEx("pfilter")
### * pfilter

flush(stderr()); flush(stdout())

### Name: pfilter
### Title: Bootstrap particle filter
### Aliases: pfilter pfilter,SimInf_model-method

### ** Examples

## Not run: 
##D ## Let us consider an SIR model in a closed population with N = 100
##D ## individuals of whom one is initially infectious and the rest are
##D ## susceptible. First, generate one realisation (with a specified
##D ## seed) from the model with known parameters 'beta = 0.16' and
##D ## 'gamma = 0.077'. Then, use 'pfilter' to apply the bootstrap
##D ## particle algorithm on the simulated data.
##D model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##D              tspan = seq(1, 100, by = 3),
##D              beta = 0.16,
##D              gamma = 0.077)
##D 
##D ## Run the SIR model to generate simulated observed data for the
##D ## number of infected individuals.
##D set.seed(22)
##D infected <- trajectory(run(model), "I")[, c("time", "I")]
##D colnames(infected) <- c("time", "Iobs")
##D 
##D ## Use a Poison observation process for the infected individuals, such
##D ## that 'Iobs ~ poison(I + 1e-6)'. A small constant '1e-6' is added to
##D ## prevent numerical errors, since the simulated counts 'I' could be
##D ## zero, which would result in the Poisson rate parameter being zero,
##D ## which violates the conditions of the Poisson distribution. Use 1000
##D ## particles.
##D pf <- pfilter(model,
##D               obs_process = Iobs ~ poisson(I + 1e-6),
##D               data = infected,
##D               n_particles = 1000)
##D 
##D ## Print a brief summary.
##D pf
##D 
##D ## Compare the number infected 'I' in the filtered trajectory with the
##D ## infected 'Iobs' in the observed data.
##D plot(pf, ~I)
##D lines(Iobs ~ time, infected, col = "blue", lwd = 2, type = "s")
## End(Not run)



cleanEx()
nameEx("plot")
### * plot

flush(stderr()); flush(stdout())

### Name: plot,SimInf_model-method
### Title: Display the outcome from a simulated trajectory
### Aliases: plot,SimInf_model-method

### ** Examples

## Not run: 
##D ## For reproducibility, call the set.seed() function and specify
##D ## the number of threads to use. To use all available threads,
##D ## remove the set_num_threads() call.
##D set.seed(123)
##D set_num_threads(1)
##D 
##D ## Create an 'SIR' model with 100 nodes and initialise
##D ## it with 990 susceptible individuals and 10 infected
##D ## individuals in each node. Run the model over 100 days.
##D model <- SIR(u0 = data.frame(S = rep(990, 100),
##D                              I = rep(10, 100),
##D                              R = rep(0, 100)),
##D              tspan = 1:100,
##D              beta = 0.16,
##D              gamma = 0.077)
##D 
##D ## Run the model and save the result.
##D result <- run(model)
##D 
##D ## Plot the median and interquartile range of the number
##D ## of susceptible, infected and recovered individuals.
##D plot(result)
##D 
##D ## Plot the median and the middle 95\##D 
##D ## number of susceptible, infected and recovered individuals.
##D plot(result, range = 0.95)
##D 
##D ## Plot the median and interquartile range of the  number
##D ## of infected individuals.
##D plot(result, "I")
##D 
##D ## Use the formula notation instead to plot the median and
##D ## interquartile range of the number of infected individuals.
##D plot(result, ~I)
##D 
##D ## Plot the number of susceptible, infected
##D ## and recovered individuals in the first
##D ## three nodes.
##D plot(result, index = 1:3, range = FALSE)
##D 
##D ## Use plot type line instead.
##D plot(result, index = 1:3, range = FALSE, type = "l")
##D 
##D ## Plot the number of infected individuals in the first node.
##D plot(result, "I", index = 1, range = FALSE)
##D 
##D ## Plot the proportion of infected individuals (cases)
##D ## in the population.
##D plot(result, I ~ S + I + R)
##D 
##D ## Plot the proportion of nodes with infected individuals.
##D plot(result, I ~ S + I + R, level = 2)
##D 
##D ## Plot the median and interquartile range of the proportion
##D ## of infected individuals in each node
##D plot(result, I ~ S + I + R, level = 3)
##D 
##D ## Plot the proportion of infected individuals in the first
##D ## three nodes.
##D plot(result, I ~ S + I + R, level = 3, index = 1:3, range = FALSE)
## End(Not run)



cleanEx()
nameEx("prevalence-SimInf_model-method")
### * prevalence-SimInf_model-method

flush(stderr()); flush(stdout())

### Name: prevalence,SimInf_model-method
### Title: Calculate prevalence from a model object with trajectory data
### Aliases: prevalence,SimInf_model-method

### ** Examples

## Create an 'SIR' model with 6 nodes and initialize
## it to run over 10 days.
u0 <- data.frame(S = 100:105, I = c(0, 1, 0, 2, 0, 3), R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Determine the proportion of infected individuals (cases)
## in the population at the time-points in 'tspan'.
prevalence(result, I ~ S + I + R)

## Identical result is obtained with the shorthand 'I~.'
prevalence(result, I ~ .)

## Determine the proportion of nodes with infected individuals at
## the time-points in 'tspan'.
prevalence(result, I ~ S + I + R, level = 2)

## Determine the proportion of infected individuals in each node
## at the time-points in 'tspan'.
prevalence(result, I ~ S + I + R, level = 3)

## Determine the proportion of infected individuals in each node
## at the time-points in 'tspan' when the number of recovered is
## zero.
prevalence(result, I ~ S + I + R | R == 0, level = 3)



cleanEx()
nameEx("punchcard-set")
### * punchcard-set

flush(stderr()); flush(stdout())

### Name: punchcard<-
### Title: Set a template for where to record result during a simulation
### Aliases: punchcard<- punchcard<-,SimInf_model-method

### ** Examples

## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 6 nodes and initialize it to run over 10 days.
u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Run the model.
result <- run(model)

## Display the trajectory with data for every node at each
## time-point in tspan.
trajectory(result)

## Assume we are only interested in nodes '2' and '4' at the
## time-points '3' and '5'
df <- data.frame(time = c(3, 5, 3, 5),
                 node = c(2, 2, 4, 4),
                 S = c(TRUE, TRUE, TRUE, TRUE),
                 I = c(TRUE, TRUE, TRUE, TRUE),
                 R = c(TRUE, TRUE, TRUE, TRUE))
punchcard(model) <- df
result <- run(model)
trajectory(result)

## We can also specify to record only some of the compartments in
## each time-step.
df <- data.frame(time = c(3, 5, 3, 5),
                 node = c(2, 2, 4, 4),
                 S = c(FALSE, TRUE, TRUE, TRUE),
                 I = c(TRUE, FALSE, TRUE, FALSE),
                 R = c(TRUE, FALSE, TRUE, TRUE))
punchcard(model) <- df
result <- run(model)
trajectory(result)

## A shortcut to specify to record all of the compartments in
## each time-step is to only inlude node and time.
df <- data.frame(time = c(3, 5, 3, 5),
                 node = c(2, 2, 4, 4))
punchcard(model) <- df
result <- run(model)
trajectory(result)

## It is possible to use an empty 'data.frame' to specify
## that no data-points should be recorded for the trajectory.
punchcard(model) <- data.frame()
result <- run(model)
trajectory(result)

## Use 'NULL' to reset the model to record data for every node at
## each time-point in tspan.
punchcard(model) <- NULL
result <- run(model)
trajectory(result)



cleanEx()
nameEx("run")
### * run

flush(stderr()); flush(stdout())

### Name: run
### Title: Run the SimInf stochastic simulation algorithm
### Aliases: run run,SimInf_model-method run,SEIR-method run,SIR-method
###   run,SIS-method run,SISe-method run,SISe3-method run,SISe3_sp-method
###   run,SISe_sp-method run,SimInf_abc-method

### ** Examples

## For reproducibility, call the set.seed() function and specify
## the number of threads to use. To use all available threads,
## remove the set_num_threads() call.
set.seed(123)
set_num_threads(1)

## Create an 'SIR' model with 10 nodes and initialise
## it to run over 100 days.
model <- SIR(u0 = data.frame(S = rep(99, 10),
                             I = rep(1, 10),
                             R = rep(0, 10)),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the model and save the result.
result <- run(model)

## Plot the proportion of susceptible, infected and recovered
## individuals.
plot(result)



cleanEx()
nameEx("select_matrix-set")
### * select_matrix-set

flush(stderr()); flush(stdout())

### Name: select_matrix<-
### Title: Set the select matrix for a 'SimInf_model' object
### Aliases: select_matrix<- select_matrix<-,SimInf_model-method

### ** Examples

## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Set the select matrix
select_matrix(model) <- matrix(c(1, 0, 0, 1, 1, 1, 0, 0, 1), nrow = 3)

## Extract the select matrix from the model
select_matrix(model)



cleanEx()
nameEx("select_matrix")
### * select_matrix

flush(stderr()); flush(stdout())

### Name: select_matrix
### Title: Extract the select matrix from a 'SimInf_model' object
### Aliases: select_matrix select_matrix,SimInf_model-method

### ** Examples

## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Extract the select matrix from the model
select_matrix(model)



cleanEx()
nameEx("shift_matrix-set")
### * shift_matrix-set

flush(stderr()); flush(stdout())

### Name: shift_matrix<-
### Title: Set the shift matrix for a 'SimInf_model' object
### Aliases: shift_matrix<- shift_matrix<-,SimInf_model-method

### ** Examples

## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Set the shift matrix
shift_matrix(model) <- matrix(c(2, 1, 0), nrow = 3)

## Extract the shift matrix from the model
shift_matrix(model)



cleanEx()
nameEx("shift_matrix")
### * shift_matrix

flush(stderr()); flush(stdout())

### Name: shift_matrix
### Title: Extract the shift matrix from a 'SimInf_model' object
### Aliases: shift_matrix shift_matrix,SimInf_model-method

### ** Examples

## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Extract the shift matrix from the model
shift_matrix(model)



cleanEx()
nameEx("show-SimInf_model-method")
### * show-SimInf_model-method

flush(stderr()); flush(stdout())

### Name: show,SimInf_model-method
### Title: Brief summary of 'SimInf_model'
### Aliases: show,SimInf_model-method

### ** Examples

## Create an 'SIR' model with 10 nodes and initialise
## it to run over 100 days.
model <- SIR(u0 = data.frame(S = rep(99, 10),
                             I = rep(1, 10),
                             R = rep(0, 10)),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Brief summary of the model
model

## Run the model and save the result
result <- run(model)

## Brief summary of the result. Note that 'U' and 'V' are
## non-empty after running the model.
result



cleanEx()
nameEx("trajectory-SimInf_model-method")
### * trajectory-SimInf_model-method

flush(stderr()); flush(stdout())

### Name: trajectory,SimInf_model-method
### Title: Extract data from a simulated trajectory
### Aliases: trajectory,SimInf_model-method

### ** Examples

## Create an 'SIR' model with 6 nodes and initialize
## it to run over 10 days.
u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)

## Run the model to generate a single stochastic trajectory.
result <- run(model)

## Extract the number of individuals in each compartment at the
## time-points in 'tspan'.
trajectory(result)

## Extract the number of recovered individuals in the first node
## at the time-points in 'tspan'.
trajectory(result, compartments = "R", index = 1)

## Extract the number of recovered individuals in the first and
## third node at the time-points in 'tspan'.
trajectory(result, compartments = "R", index = c(1, 3))

## Create an 'SISe' model with 6 nodes and initialize
## it to run over 10 days.
u0 <- data.frame(S = 100:105, I = 1:6)
model <- SISe(u0 = u0, tspan = 1:10, phi = rep(0, 6),
    upsilon = 0.02, gamma = 0.1, alpha = 1, epsilon = 1.1e-5,
    beta_t1 = 0.15, beta_t2 = 0.15, beta_t3 = 0.15, beta_t4 = 0.15,
    end_t1 = 91, end_t2 = 182, end_t3 = 273, end_t4 = 365)

## Run the model
result <- run(model)

## Extract the continuous state variable 'phi' which represents
## the environmental infectious pressure.
trajectory(result, "phi")



cleanEx()
nameEx("u0-set")
### * u0-set

flush(stderr()); flush(stdout())

### Name: u0<-
### Title: Update the initial compartment state u0 in each node
### Aliases: u0<- u0<-,SimInf_model-method

### ** Examples

## Create an SIR model object.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the SIR model and plot the result.
set.seed(22)
result <- run(model)
plot(result)

## Update u0 and run the model again
u0(model) <- data.frame(S = 990, I = 10, R = 0)
result <- run(model)
plot(result)



cleanEx()
nameEx("u0")
### * u0

flush(stderr()); flush(stdout())

### Name: u0
### Title: Get the initial compartment state
### Aliases: u0 u0,SimInf_model-method u0,SimInf_indiv_events-method

### ** Examples

## Create an SIR model object.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Get the initial compartment state.
u0(model)



cleanEx()
nameEx("u0_SEIR")
### * u0_SEIR

flush(stderr()); flush(stdout())

### Name: u0_SEIR
### Title: Example data to initialize the 'SEIR' model
### Aliases: u0_SEIR

### ** Examples

## Not run: 
##D ## For reproducibility, call the set.seed() function and specify
##D ## the number of threads to use. To use all available threads,
##D ## remove the set_num_threads() call.
##D set.seed(123)
##D set_num_threads(1)
##D 
##D ## Create an 'SEIR' model with 1600 nodes and initialize it to
##D ## run over 4*365 days and record data at weekly time-points.
##D ## Add ten infected individuals to the first node.
##D u0 <- u0_SEIR()
##D u0$I[1] <- 10
##D tspan <- seq(from = 1, to = 4*365, by = 7)
##D model <- SEIR(u0      = u0,
##D               tspan   = tspan,
##D               events  = events_SEIR(),
##D               beta    = 0.16,
##D               epsilon = 0.25,
##D               gamma   = 0.01)
##D 
##D ## Run the model to generate a single stochastic trajectory.
##D result <- run(model)
##D plot(result)
##D 
##D ## Summarize trajectory
##D summary(result)
## End(Not run)



cleanEx()
nameEx("u0_SIR")
### * u0_SIR

flush(stderr()); flush(stdout())

### Name: u0_SIR
### Title: Example data to initialize the 'SIR' model
### Aliases: u0_SIR

### ** Examples

## Not run: 
##D ## For reproducibility, call the set.seed() function and specify
##D ## the number of threads to use. To use all available threads,
##D ## remove the set_num_threads() call.
##D set.seed(123)
##D set_num_threads(1)
##D 
##D ## Create an 'SIR' model with 1600 nodes and initialize
##D ## it to run over 4*365 days. Add one infected individual
##D ## to the first node.
##D u0 <- u0_SIR()
##D u0$I[1] <- 1
##D tspan <- seq(from = 1, to = 4*365, by = 1)
##D model <- SIR(u0     = u0,
##D              tspan  = tspan,
##D              events = events_SIR(),
##D              beta   = 0.16,
##D              gamma  = 0.01)
##D 
##D ## Run the model to generate a single stochastic trajectory.
##D result <- run(model)
##D plot(result)
##D 
##D ## Summarize trajectory
##D summary(result)
## End(Not run)



cleanEx()
nameEx("u0_SIS")
### * u0_SIS

flush(stderr()); flush(stdout())

### Name: u0_SIS
### Title: Example data to initialize the 'SIS' model
### Aliases: u0_SIS

### ** Examples

## Not run: 
##D ## For reproducibility, call the set.seed() function and specify
##D ## the number of threads to use. To use all available threads,
##D ## remove the set_num_threads() call.
##D set.seed(123)
##D set_num_threads(1)
##D 
##D ## Create an 'SIS' model with 1600 nodes and initialize
##D ## it to run over 4*365 days. Add one infected individual
##D ## to the first node.
##D u0 <- u0_SIS()
##D u0$I[1] <- 1
##D tspan <- seq(from = 1, to = 4*365, by = 1)
##D model <- SIS(u0     = u0,
##D              tspan  = tspan,
##D              events = events_SIS(),
##D              beta   = 0.16,
##D              gamma  = 0.01)
##D 
##D ## Run the model to generate a single stochastic trajectory.
##D result <- run(model)
##D plot(result)
##D 
##D ## Summarize trajectory
##D summary(result)
## End(Not run)



cleanEx()
nameEx("u0_SISe")
### * u0_SISe

flush(stderr()); flush(stdout())

### Name: u0_SISe
### Title: Example data to initialize the 'SISe' model
### Aliases: u0_SISe

### ** Examples

## Not run: 
##D ## For reproducibility, call the set.seed() function and specify
##D ## the number of threads to use. To use all available threads,
##D ## remove the set_num_threads() call.
##D set.seed(123)
##D set_num_threads(1)
##D 
##D ## Create an 'SISe' model with 1600 nodes and initialize it to
##D ## run over 4*365 days and record data at weekly time-points.
##D 
##D ## Load the initial population and add ten infected individuals to
##D ## the first node.
##D u0 <- u0_SISe()
##D u0$I[1] <- 10
##D 
##D ## Define 'tspan' to run the simulation over 4*365 and record the
##D ## state of the system at weekly time-points.
##D tspan <- seq(from = 1, to = 4*365, by = 7)
##D 
##D ## Load scheduled events for the population of nodes with births,
##D ## deaths and between-node movements of individuals.
##D events <- events_SISe()
##D 
##D ## Create an 'SISe' model
##D model <- SISe(u0 = u0, tspan = tspan, events = events_SISe(),
##D               phi = 0, upsilon = 1.8e-2, gamma = 0.1, alpha = 1,
##D               beta_t1 = 1.0e-1, beta_t2 = 1.0e-1, beta_t3 = 1.25e-1,
##D               beta_t4 = 1.25e-1, end_t1 = 91, end_t2 = 182,
##D               end_t3 = 273, end_t4 = 365, epsilon = 0)
##D 
##D ## Run the model to generate a single stochastic trajectory.
##D result <- run(model)
##D 
##D ## Summarize trajectory
##D summary(result)
##D 
##D ## Plot the proportion of nodes with at least one infected
##D ## individual.
##D plot(result, I~S+I, level = 2, type = "l")
## End(Not run)



cleanEx()
nameEx("u0_SISe3")
### * u0_SISe3

flush(stderr()); flush(stdout())

### Name: u0_SISe3
### Title: Example data to initialize the 'SISe3' model
### Aliases: u0_SISe3
### Keywords: dataset

### ** Examples

## Not run: 
##D ## For reproducibility, call the set.seed() function and specify
##D ## the number of threads to use. To use all available threads,
##D ## remove the set_num_threads() call.
##D set.seed(123)
##D set_num_threads(1)
##D 
##D ## Create an 'SISe3' model with 1600 nodes and initialize it to
##D ## run over 4*365 days and record data at weekly time-points.
##D 
##D ## Load the initial population and add ten infected individuals to
##D ## I_1 in the first node.
##D u0 <- u0_SISe3
##D u0$I_1[1] <- 10
##D 
##D ## Define 'tspan' to run the simulation over 4*365 and record the
##D ## state of the system at weekly time-points.
##D tspan <- seq(from = 1, to = 4*365, by = 7)
##D 
##D ## Load scheduled events for the population of nodes with births,
##D ## deaths and between-node movements of individuals.
##D events <- events_SISe3
##D 
##D ## Create a 'SISe3' model
##D model <- SISe3(u0 = u0, tspan = tspan, events = events,
##D                phi = rep(0, nrow(u0)), upsilon_1 = 1.8e-2,
##D                upsilon_2 = 1.8e-2, upsilon_3 = 1.8e-2,
##D                gamma_1 = 0.1, gamma_2 = 0.1, gamma_3 = 0.1,
##D                alpha = 1, beta_t1 = 1.0e-1, beta_t2 = 1.0e-1,
##D                beta_t3 = 1.25e-1, beta_t4 = 1.25e-1, end_t1 = 91,
##D                end_t2 = 182, end_t3 = 273, end_t4 = 365, epsilon = 0)
##D 
##D ## Run the model to generate a single stochastic trajectory.
##D result <- run(model)
##D 
##D ## Summarize trajectory
##D summary(result)
##D 
##D ## Plot the proportion of nodes with at least one infected
##D ## individual.
##D plot(result, I_1 + I_2 + I_3 ~ ., level = 2, type = "l")
## End(Not run)



cleanEx()
nameEx("v0-set")
### * v0-set

flush(stderr()); flush(stdout())

### Name: v0<-
### Title: Update the initial continuous state v0 in each node
### Aliases: v0<- v0<-,SimInf_model-method

### ** Examples

## Create an 'SISe' model with no infected individuals and no
## infectious pressure (phi = 0, epsilon = 0).
model <- SISe(u0 = data.frame(S = 100, I = 0), tspan = 1:100,
              phi = 0, upsilon = 0.02, gamma = 0.1, alpha = 1,
              epsilon = 0, beta_t1 = 0.15, beta_t2 = 0.15,
              beta_t3 = 0.15, beta_t4 = 0.15, end_t1 = 91,
              end_t2 = 182, end_t3 = 273, end_t4 = 365)

## Run the 'SISe' model and plot the result.
set.seed(22)
result <- run(model)
plot(result)

## Update the infectious pressure 'phi' in 'v0' and run
## the model again.
v0(model) <- data.frame(phi = 1)
result <- run(model)
plot(result)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
