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

##' Definition of the \code{SISe} model
##'
##' Class to handle the SISe \code{\link{SimInf_model}}.
##' @include SimInf_model.R
##' @export
setClass("SISe", contains = c("SimInf_model"))

##' Create a SISe model
##'
##' Create an \sQuote{SISe} model to be used by the simulation
##' framework.
##'
##' The \sQuote{SISe} model contains two compartments; number of
##' susceptible (S) and number of infectious (I). Additionally, it
##' contains an environmental compartment to model shedding of a
##' pathogen to the environment. Consequently, the model has two state
##' transitions,
##'
##' \deqn{S \stackrel{\upsilon \varphi S}{\longrightarrow} I}{
##' S -- upsilon phi S --> I}
##'
##' \deqn{I \stackrel{\gamma I}{\longrightarrow} S}{
##' I -- gamma I --> S}
##'
##' where the transition rate per unit of time from susceptible to
##' infected is proportional to the concentration of the environmental
##' contamination \eqn{\varphi}{phi} in each node. Moreover, the
##' transition rate from infected to susceptible is the recovery rate
##' \eqn{\gamma}, measured per individual and per unit of
##' time. Finally, the environmental infectious pressure in each node
##' is evolved by,
##'
##' \deqn{\frac{d\varphi(t)}{dt} = \frac{\alpha I(t)}{N(t)} - \beta(t)
##' \varphi(t) + \epsilon}{
##' dphi(t) / dt = alpha I(t) / N(t) - beta(t) phi(t) + epsilon}
##'
##' where \eqn{\alpha} is the average shedding rate of the pathogen to
##' the environment per infected individual and \eqn{N = S + I} the
##' size of the node. The seasonal decay and removal of the pathogen
##' is captured by \eqn{\beta(t)}. It is also possible to include a
##' small background infectious pressure \eqn{\epsilon} to allow for
##' other indirect sources of environmental contamination. The
##' environmental infectious pressure \eqn{\varphi(t)}{phi(t)} in each
##' node is evolved each time unit by the Euler forward method. The
##' value of \eqn{\varphi(t)}{phi(t)} is saved at the time-points
##' specified in \code{tspan}.
##'
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of sucsceptible in each node}
##' \item{I}{The number of infected in each node}
##' }
##'
##' @template beta-section
##' @template u0-param
##' @template tspan-param
##' @template events-param
##' @template phi-param
##' @param upsilon Indirect transmission rate of the environmental
##'     infectious pressure
##' @param gamma The recovery rate from infected to susceptible
##' @param alpha Shed rate from infected individuals
##' @template beta-end-param
##' @param epsilon The background environmental infectious pressure
##' @return \code{SISe}
##' @include check_arguments.R
##' @export
##' @importFrom methods as
SISe <- function(u0,
                 tspan,
                 events  = NULL,
                 phi     = NULL,
                 upsilon = NULL,
                 gamma   = NULL,
                 alpha   = NULL,
                 beta_t1 = NULL,
                 beta_t2 = NULL,
                 beta_t3 = NULL,
                 beta_t4 = NULL,
                 end_t1  = NULL,
                 end_t2  = NULL,
                 end_t3  = NULL,
                 end_t4  = NULL,
                 epsilon = NULL) {
    compartments <- c("S", "I")

    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments)

    ## Check initial infectious pressure
    if (is.null(phi))
        phi <- 0
    phi <- rep(phi, length.out = nrow(u0))
    check_infectious_pressure_arg(nrow(u0), phi)

    ## Check for non-numeric parameters
    check_gdata_arg(upsilon, gamma, alpha, beta_t1, beta_t2, beta_t3, beta_t4,
                    epsilon)

    ## Check interval endpoints
    check_integer_arg(end_t1, end_t2, end_t3, end_t4)
    end_t1 <- rep(end_t1, length.out = nrow(u0))
    end_t2 <- rep(end_t2, length.out = nrow(u0))
    end_t3 <- rep(end_t3, length.out = nrow(u0))
    end_t4 <- rep(end_t4, length.out = nrow(u0))
    check_end_t_arg(nrow(u0), end_t1, end_t2, end_t3, end_t4)

    ## Arguments seem ok...go on

    E <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2,
                dimnames = list(compartments, c("1", "2")))

    G <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2,
                dimnames = list(c("S -> upsilon*phi*S -> I",
                                  "I -> gamma*I -> S"),
                                c("1", "2")))

    S <- matrix(c(-1,  1, 1, -1), nrow = 2, ncol = 2,
                dimnames = list(compartments, c("1", "2")))

    v0 <- matrix(as.numeric(phi), nrow  = 1, byrow = TRUE,
                 dimnames = list("phi"))

    ldata <- matrix(as.numeric(c(end_t1, end_t2, end_t3, end_t4)),
                    nrow  = 4, byrow = TRUE,
                    dimnames = list(c("end_t1", "end_t2", "end_t3", "end_t4")))

    gdata <- as.numeric(c(upsilon, gamma, alpha, beta_t1, beta_t2,
                          beta_t3, beta_t4, epsilon))
    names(gdata) <- c("upsilon", "gamma", "alpha", "beta_t1", "beta_t2",
                      "beta_t3", "beta_t4", "epsilon")

    model <- SimInf_model(G      = G,
                          S      = S,
                          E      = E,
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0)

    as(model, "SISe")
}

##' Example data to initialize events for the \sQuote{SISe} model
##'
##' Example data to initialize scheduled events for a population of
##' 1600 nodes and demonstrate the \code{\linkS4class{SISe}} model.
##'
##' Example data to initialize scheduled events (see
##' \code{\linkS4class{SimInf_events}}) for a population of 1600 nodes
##' and demonstrate the \code{\linkS4class{SISe}} model. The dataset
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
##' ## Create an 'SISe' model with 1600 nodes and initialize
##' ## it to run over 4*365 days. Add one infected individual
##' ## to the first node.
##' u0 <- u0_SISe()
##' u0$I[1] <- 1
##' tspan <- seq(from = 1, to = 4*365, by = 1)
##' model <- SISe(u0 = u0, tspan = tspan, events = events_SISe(),
##'               phi = 0, upsilon = 1.8e-2, gamma = 0.1, alpha = 1,
##'               beta_t1 = 1.0e-1, beta_t2 = 1.0e-1, beta_t3 = 1.25e-1,
##'               beta_t4 = 1.25e-1, end_t1 = 91, end_t2 = 182,
##'               end_t3 = 273, end_t4 = 365, epsilon = 0)
##'
##' ## Display the number of individuals affected by each event type
##' ## per day.
##' plot(events(model))
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##'
##' ## Summarize the trajectory. The summary includes the number of
##' ## events by event type.
##' summary(result)
events_SISe <- function() {
    data("events_SISe3", package = "SimInf", envir = environment())
    events_SISe3$select[events_SISe3$event == "exit"] <- 2L
    events_SISe3$select[events_SISe3$event == "enter"] <- 1L
    events_SISe3 <- events_SISe3[events_SISe3$event != "intTrans", ]
    events_SISe3$select[events_SISe3$event == "extTrans"] <- 2L
    events_SISe3
}

##' Example data to initialize the \sQuote{SISe} model
##'
##' Example data to initialize a population of 1600 nodes and
##' demonstrate the \code{\linkS4class{SISe}} model.
##'
##' A \code{data.frame} with the number of individuals in the
##' \sQuote{S} and \sQuote{I} compartments in 1600 nodes. Note that
##' the \sQuote{I} compartment is zero.
##' @return A \code{data.frame}
##' @export
##' @importFrom utils data
##' @examples
##' ## Create an 'SISe' model with 1600 nodes and initialize it to
##' ## run over 4*365 days and record data at weekly time-points.
##'
##' ## Load the initial population and add ten infected individuals to
##' ## the first node.
##' u0 <- u0_SISe()
##' u0$I[1] <- 10
##'
##' ## Define 'tspan' to run the simulation over 4*365 and record the
##' ## state of the system at weekly time-points.
##' tspan <- seq(from = 1, to = 4*365, by = 7)
##'
##' ## Load scheduled events for the population of nodes with births,
##' ## deaths and between-node movements of individuals.
##' events <- events_SISe()
##'
##' ## Create an 'SISe' model
##' model <- SISe(u0 = u0, tspan = tspan, events = events_SISe(),
##'               phi = 0, upsilon = 1.8e-2, gamma = 0.1, alpha = 1,
##'               beta_t1 = 1.0e-1, beta_t2 = 1.0e-1, beta_t3 = 1.25e-1,
##'               beta_t4 = 1.25e-1, end_t1 = 91, end_t2 = 182,
##'               end_t3 = 273, end_t4 = 365, epsilon = 0)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##'
##' ## Summarize trajectory
##' summary(result)
##'
##' ## Plot the proportion of nodes with at least one infected
##' ## individual.
##' plot(result, I~S+I, level = 2, type = "l")
u0_SISe <- function() {
    u0 <- u0_SIR()
    u0[, c("S", "I")]
}
