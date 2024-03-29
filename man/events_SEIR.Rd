% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEIR.R
\name{events_SEIR}
\alias{events_SEIR}
\title{Example data to initialize events for the \sQuote{SEIR} model}
\usage{
events_SEIR()
}
\value{
A \code{data.frame}
}
\description{
Example data to initialize scheduled events for a population of
1600 nodes and demonstrate the \code{\linkS4class{SEIR}} model.
}
\details{
Example data to initialize scheduled events (see
\code{\linkS4class{SimInf_events}}) for a population of 1600 nodes
and demonstrate the \code{\linkS4class{SEIR}} model. The dataset
contains 466692 events for 1600 nodes distributed over 4 * 365
days. The events are divided into three types: \sQuote{Exit}
events remove individuals from the population (n = 182535),
\sQuote{Enter} events add individuals to the population (n =
182685), and \sQuote{External transfer} events move individuals
between nodes in the population (n = 101472). The vignette
contains a detailed description of how scheduled events operate on
a model.
}
\examples{
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
}
