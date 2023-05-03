## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2023 Stefan Widgren
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

##' Definition of the \acronym{SIS} model
##'
##' Class to handle the \acronym{SIS} \code{\link{SimInf_model}}.
##'
##' The \acronym{SIS} model contains two compartments; number of
##' susceptible (S), and number of infectious (I).  Moreover, it has
##' two state transitions, \deqn{S \stackrel{\beta S I /
##' N}{\longrightarrow} I}{ S -- beta S I / N --> I} \deqn{I
##' \stackrel{\gamma I}{\longrightarrow} S}{I -- gamma I --> S} where
##' \eqn{\beta} is the transmission rate, \eqn{\gamma} is the recovery
##' rate, and \eqn{N=S+I}.
##' @include SimInf_model.R
##' @export
##' @examples
##' ## Create an SIS model object.
##' model <- SIS(u0 = data.frame(S = 99, I = 1),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the SIS model and plot the result.
##' set.seed(22)
##' result <- run(model)
##' plot(result)
setClass("SIS", contains = c("SimInf_model"))

##' The compartments in an SIS model
##' @noRd
compartments_SIS <- function()
    c("S", "I")
}

##' The select matrix 'E' for an SIS model
##' @noRd
select_matrix_SIS <- function() {
    matrix(c(1, 0, 1, 1),
           nrow = 2,
           ncol = 2,
           dimnames = list(compartments_SIS(), seq_len(2)))
}

##' Create an \acronym{SIS} model
##'
##' Create an \acronym{SIS} model to be used by the simulation
##' framework.
##'
##' The \acronym{SIS} model contains two compartments; number of
##' susceptible (S), and number of infectious (I).  Moreover, it has
##' two state transitions, \deqn{S \stackrel{\beta S I /
##' N}{\longrightarrow} I}{ S -- beta S I / N --> I} \deqn{I
##' \stackrel{\gamma I}{\longrightarrow} S}{I -- gamma I --> S} where
##' \eqn{\beta} is the transmission rate, \eqn{\gamma} is the recovery
##' rate, and \eqn{N=S+I}.
##'
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of sucsceptible in each node}
##' \item{I}{The number of infected in each node}
##' }
##'
##' @template u0-param
##' @template tspan-param
##' @template events-param
##' @template beta-param
##' @template gamma-param
##' @return A \code{\link{SimInf_model}} of class \code{SIS}
##' @include check_arguments.R
##' @export
##' @examples
##' ## Create an SIS model object.
##' model <- SIS(u0 = data.frame(S = 99, I = 1),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the SIS model and plot the result.
##' set.seed(22)
##' result <- run(model)
##' plot(result)
SIS <- function(u0,
                tspan,
                events = NULL,
                beta   = NULL,
                gamma  = NULL) {
    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments_SIS())

    ## Check for non-numeric parameters
    check_ldata_arg(nrow(u0), beta, gamma)
    beta <- rep(beta, length.out = nrow(u0))
    gamma <- rep(gamma, length.out = nrow(u0))

    ## Arguments seem ok...go on

    G <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2,
                dimnames = list(c("S -> upsilon*S*I -> I",
                                  "I -> gamma*I -> S"),
                                c("1", "2")))

    S <- matrix(c(-1,  1, 1, -1), nrow = 2, ncol = 2,
                dimnames = list(compartments_SIS(), c("1", "2")))

    ldata <- matrix(as.numeric(c(beta, gamma)),
                    nrow  = 2, byrow = TRUE,
                    dimnames = list(c("beta", "gamma")))

    model <- SimInf_model(G      = G,
                          S      = S,
                          E      = select_matrix_SIS(),
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          u0     = u0)

    methods::as(model, "SIS")
}

##' Example data to initialize events for the \sQuote{SIS} model
##'
##' Example data to initialize scheduled events for a population of
##' 1600 nodes and demonstrate the \code{\linkS4class{SIS}} model.
##'
##' Example data to initialize scheduled events (see
##' \code{\linkS4class{SimInf_events}}) for a population of 1600 nodes
##' and demonstrate the \code{\linkS4class{SIS}} model. The dataset
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
##' @examples
##' ## Create an 'SIS' model with 1600 nodes and initialize
##' ## it to run over 4*365 days. Add one infected individual
##' ## to the first node.
##' u0 <- u0_SIS()
##' u0$I[1] <- 1
##' tspan <- seq(from = 1, to = 4*365, by = 1)
##' model <- SIS(u0     = u0,
##'              tspan  = tspan,
##'              events = events_SIS(),
##'              beta   = 0.16,
##'              gamma  = 0.01)
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
events_SIS <- function() {
    events_SISe()
}

##' Example data to initialize the \sQuote{SIS} model
##'
##' Example data to initialize a population of 1600 nodes and
##' demonstrate the \code{\linkS4class{SIS}} model.
##'
##' A \code{data.frame} with the number of individuals in the
##' \sQuote{S}, and \sQuote{I} compartments in 1600 nodes. Note that
##' the \sQuote{I} compartment is zero.
##' @return A \code{data.frame}
##' @export
##' @examples
##' ## Create an 'SIS' model with 1600 nodes and initialize
##' ## it to run over 4*365 days. Add one infected individual
##' ## to the first node.
##' u0 <- u0_SIS()
##' u0$I[1] <- 1
##' tspan <- seq(from = 1, to = 4*365, by = 1)
##' model <- SIS(u0     = u0,
##'              tspan  = tspan,
##'              events = events_SIS(),
##'              beta   = 0.16,
##'              gamma  = 0.01)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##' plot(result)
##'
##' ## Summarize trajectory
##' summary(result)
u0_SIS <- function() {
    u0_SISe()
}
