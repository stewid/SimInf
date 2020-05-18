## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2020 Stefan Widgren
##
## SimInf is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SimInf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

##' Definition of the \sQuote{SEIR} model
##'
##' Class to handle the SEIR \code{\link{SimInf_model}}.
##' @include SimInf_model.R
##' @export
setClass("SEIR", contains = c("SimInf_model"))

##' Create an \acronym{SEIR} model
##'
##' Create an \acronym{SEIR} model to be used by the simulation
##' framework.
##'
##' The \acronym{SEIR} model contains four compartments; number of
##' susceptible (S), number of exposed (E) (those who have been
##' infected but are not yet infectious), number of infectious (I),
##' and number of recovered (R).  Moreover, it has three state
##' transitions,
##'
##' \deqn{S \stackrel{\beta S I / N}{\longrightarrow} E}{
##'   S -- beta S I / N --> E}
##' \deqn{E \stackrel{\epsilon E}{\longrightarrow} I}{E -- epsilon E --> I}
##' \deqn{I \stackrel{\gamma I}{\longrightarrow} R}{I -- gamma I --> R}
##'
##' where \eqn{\beta} is the transmission rate, \eqn{\epsilon} is the
##' incubation rate, \eqn{\gamma} is the recovery rate, and
##' \eqn{N=S+E+I+R}.
##'
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of sucsceptible in each node}
##' \item{E}{The number of exposed in each node}
##' \item{I}{The number of infected in each node}
##' \item{R}{The number of recovered in each node}
##' }
##'
##' @template u0-param
##' @template tspan-param
##' @template events-param
##' @template beta-param
##' @template epsilon-param
##' @template gamma-param
##' @return A \code{\link{SimInf_model}} of class \code{SEIR}
##' @include check_arguments.R
##' @export
##' @examples
##' ## Create a SEIR model object.
##' model <- SEIR(u0 = data.frame(S = 99, E = 0, I = 1, R = 0),
##'               tspan = 1:100,
##'               beta = 0.16,
##'               epsilon = 0.25,
##'               gamma = 0.077)
##'
##' ## Run the SEIR model and plot the result.
##' set.seed(3)
##' result <- run(model)
##' plot(result)
SEIR <- function(u0,
                 tspan,
                 events  = NULL,
                 beta    = NULL,
                 epsilon = NULL,
                 gamma   = NULL) {
    compartments <- c("S", "E", "I", "R")

    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments)

    ## Check for non-numeric parameters
    check_ldata_arg(nrow(u0), beta, epsilon, gamma)
    beta <- rep(beta, length.out = nrow(u0))
    epsilon <- rep(epsilon, length.out = nrow(u0))
    gamma <- rep(gamma, length.out = nrow(u0))

    ## Arguments seem ok...go on

    E <- matrix(c(1, 0, 0, 0, 1, 1, 1, 1), nrow = 4, ncol = 2,
                dimnames = list(compartments, c("1", "2")))

    G <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1), nrow = 3, ncol = 3,
                dimnames = list(c("S -> beta*S*I/(S+E+I+R) -> E",
                                  "E -> epsilon*E -> I",
                                  "I -> gamma*I -> R"),
                                c("1", "2", "3")))

    S <- matrix(c(-1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1), nrow = 4, ncol = 3,
                dimnames = list(compartments, c("1", "2", "3")))

    ldata <- matrix(as.numeric(c(beta, epsilon, gamma)),
                    nrow  = 3, byrow = TRUE,
                    dimnames = list(c("beta", "epsilon", "gamma")))

    model <- SimInf_model(G      = G,
                          S      = S,
                          E      = E,
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          u0     = u0)

    as(model, "SEIR")
}

##' Example data to initialize events for the \sQuote{SEIR} model
##'
##' Example data to initialize scheduled events for a population of
##' 1600 nodes and demonstrate the \code{\linkS4class{SEIR}} model.
##'
##' Example data to initialize scheduled events (see
##' \code{\linkS4class{SimInf_events}}) for a population of 1600 nodes
##' and demonstrate the \code{\linkS4class{SEIR}} model. The dataset
##' contains 466692 events for 1600 nodes distributed over 4 * 365
##' days. The events are divided into three types: \sQuote{Exit}
##' events remove individuals from the population (n = 182535),
##' \sQuote{Enter} events add individuals to the population (n =
##' 182685), and \sQuote{External transfer} events move individuals
##' between nodes in the population (n = 101472). The vignette
##' contains a detailed description of how scheduled events operate on
##' a model.
##' @return A \code{data.frame}
##' @export
##' @importFrom utils data
##' @examples
##' ## Create an 'SEIR' model with 1600 nodes and initialize
##' ## it to run over 4*365 days. Add one infected individual
##' ## to the first node.
##' u0 <- u0_SEIR()
##' u0$I[1] <- 1
##' tspan <- seq(from = 1, to = 4*365, by = 1)
##' model <- SEIR(u0      = u0,
##'               tspan   = tspan,
##'               events  = events_SEIR(),
##'               beta    = 0.16,
##'               epsilon = 0.25,
##'               gamma   = 0.01)
##'
##' ## Display the number of individuals affected by each event type
##' ## per day.
##' plot(events(model))
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##' plot(result)
##'
##' ## Summarize the trajectory. The summary includes the number of
##' ## events by event type.
##' summary(result)
events_SEIR <- function() {
    data("events_SISe3", package = "SimInf", envir = environment())
    events_SISe3$select[events_SISe3$event == "exit"] <- 2
    events_SISe3$select[events_SISe3$event == "enter"] <- 1
    events_SISe3 <- events_SISe3[events_SISe3$event != "intTrans", ]
    events_SISe3$select[events_SISe3$event == "extTrans"] <- 2
    events_SISe3
}

##' Example data to initialize the \sQuote{SEIR} model
##'
##' Example data to initialize a population of 1600 nodes and
##' demonstrate the \code{\linkS4class{SEIR}} model.
##'
##' A \code{data.frame} with the number of individuals in the
##' \sQuote{S}, \sQuote{E}, \sQuote{I} and \sQuote{R} compartments in
##' 1600 nodes. Note that the \sQuote{E}, \sQuote{I} and \sQuote{R}
##' compartments are zero.
##' @return A \code{data.frame}
##' @export
##' @importFrom utils data
##' @examples
##' ## Create an 'SEIR' model with 1600 nodes and initialize it to
##' ## run over 4*365 days and record data at weekly time-points.
##' ## Add ten infected individuals to the first node.
##' u0 <- u0_SEIR()
##' u0$I[1] <- 10
##' tspan <- seq(from = 1, to = 4*365, by = 7)
##' model <- SEIR(u0      = u0,
##'               tspan   = tspan,
##'               events  = events_SEIR(),
##'               beta    = 0.16,
##'               epsilon = 0.25,
##'               gamma   = 0.01)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##' plot(result)
##'
##' ## Summarize trajectory
##' summary(result)
u0_SEIR <- function() {
    u0 <- u0_SIR()
    u0$E <- 0
    u0[, c("S", "E", "I", "R")]
}
