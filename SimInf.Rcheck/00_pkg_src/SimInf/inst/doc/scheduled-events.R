## ----pressure, echo=FALSE, fig.align="left", fig.cap="**Figure 1.** Illustration of movements between nodes. Each time step depicts movements during one time unit, for example, a day. The network has *N=4* nodes where node *1* is infected and nodes *2*--*4* are non-infected. Arrows indicate movements of individuals from a source node to a destination node and labels denote the size of the shipment. Here, infection may spread from node *1* to node *3* at *t=2* and then from node *3* to node *2* at *t=3*.", out.width = '100%'----
knitr::include_graphics("img/temporal-network.svg")

## ----eval = TRUE, echo = TRUE, message = FALSE----------------------
events <- data.frame(
  event      = rep("extTrans", 6),  ## Event "extTrans" is
                                    ##  a movement between nodes
  time       = c(1, 1, 2, 2, 3, 3), ## The time that the event happens
  node       = c(3, 3, 1, 4, 3, 4), ## In which node does the event occur
  dest       = c(4, 2, 3, 3, 2, 2), ## Which node is the destination node
  n          = c(9, 2, 8, 3, 5, 4), ## How many individuals are moved
  proportion = c(0, 0, 0, 0, 0, 0), ## This is not used when n > 0
  select     = c(4, 4, 4, 4, 4, 4), ## Use the 4th column in
                                    ## the model select matrix
  shift      = c(0, 0, 0, 0, 0, 0)) ## Not used in this example

## -------------------------------------------------------------------
events

## -------------------------------------------------------------------
library(SimInf)

model <- SIR(u0 = data.frame(S = c(10, 15, 20, 25),
                             I = c(5,  0,  0,  0),
                             R = c(0,  0,  0,  0)),
             tspan = 0:3,
             beta = 0,
             gamma = 0,
             events = events)

## -------------------------------------------------------------------
model@events@E

## -------------------------------------------------------------------
set.seed(1)
set_num_threads(1)
result <- run(model)

## ----fig.width=7, fig.height=4, fig.align="left", fig.cap="**Figure 2.** Number of susceptible, infected and recovered individuals in each node."----
plot(result, range = FALSE)

## -------------------------------------------------------------------
trajectory(result)

## -------------------------------------------------------------------
u0 <- data.frame(S = c(100, 0),
                 I = c(100, 0),
                 R = c(100, 0))

## -------------------------------------------------------------------
events <- data.frame(
  event = rep("extTrans", 300), ## Event "extTrans" is a movement between nodes
  time = 1:300,                 ## The time that the event happens
  node = 1,                     ## In which node does the event occur
  dest = 2,                     ## Which node is the destination node
  n = 1,                        ## How many individuals are moved
  proportion = 0,               ## This is not used when n > 0
  select = 4,                   ## Use the 4th column in the model select matrix
  shift = 0)                    ## Not used in this example

## -------------------------------------------------------------------
model <- SIR(u0 = u0,
             tspan = 1:300,
             events = events,
             beta = 0,
             gamma = 0)

## ----fig.width=7, fig.height=4, fig.align="left", fig.cap="**Figure 3.** The individuals have an equal probability of being selected regardless of compartment."----
plot(run(model), index = 2)

## ----fig.width=7, fig.height=4, fig.align="left", fig.cap="**Figure 4.** The individuals in the $I$ compartment are more likely of being selected for a movement event."----
model@events@E[2, 4] <- 2
plot(run(model), index = 2)

## ----fig.width=7, fig.height=4, fig.align="left", fig.cap="**Figure 5.** The individuals in the $I$ compartment are even more likely of being selected for a movement event compared to the previous example."----
model@events@E[2, 4] <- 10
plot(run(model), index = 2)

## ----fig.width=7, fig.height=4, fig.align="left", fig.cap="**Figure 6.** The individuals in the $I$ and $R$ compartments are more likely of being selected for a movement event compared to individuals in the $S$ compartment."----
model@events@E[3, 4] <- 4
plot(run(model), index = 2)

