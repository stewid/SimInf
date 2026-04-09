## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2026 Stefan Widgren
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

##' Class SEIR
##'
##' Class to handle the \acronym{SEIR} model. This class inherits from
##' \code{\linkS4class{SimInf_model}}, meaning that \acronym{SEIR}
##' objects are fully compatible with all generic functions defined
##' for \code{SimInf_model}, such as \code{\link{run}},
##' \code{\link{plot}}, \code{\link{trajectory}}, and
##' \code{\link{prevalence}}.
##'
##' @template SEIR-details
##'
##' @seealso
##' \code{\link{SEIR}} for creating an \acronym{SEIR} model object,
##' \code{\linkS4class{SimInf_model}} for the parent class definition,
##' and \code{\link{SIR}} for the base model without the latent
##' period.
##' @include SimInf_model.R
##' @export
setClass("SEIR", contains = c("SimInf_model"))

##' The compartments in an SEIR model
##' @noRd
compartments_SEIR <- function() {
    c("S", "E", "I", "R")
}

##' Select matrix for events in the \acronym{SEIR} model
##'
##' Internal function returning the 4x2 select matrix (E) that maps
##' SEIR compartments (rows) to event types (columns) for event
##' processing.
##'
##' @return A 4x2 numeric matrix with compartments as rows and event
##'     types as columns. Used internally by SimInf_events.
##' @noRd
select_matrix_SEIR <- function() {
    matrix(c(1, 0, 0, 0, 1, 1, 1, 1),
           nrow = 4,
           ncol = 2,
           dimnames = list(compartments_SEIR(), seq_len(2)))
}

##' Create an \acronym{SEIR} model
##'
##' Create an \acronym{SEIR} model to be used by the simulation
##' framework.
##'
##' @template SEIR-details
##' @details
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of susceptible individuals in each node}
##' \item{E}{The number of exposed individuals in each node}
##' \item{I}{The number of infected individuals in each node}
##' \item{R}{The number of recovered individuals in each node}
##' }
##'
##' @template u0-param
##' @template tspan-param
##' @template events-param
##' @template beta-param
##' @param epsilon A numeric vector with the incubation rate from
##'     exposed to infected. Each node can have a different value. The
##'     vector must have length 1 or \code{nrow(u0)}. If the vector
##'     has length 1 but the model contains more nodes, the value is
##'     repeated for all nodes.
##' @template gamma-param
##' @return A \code{\link{SimInf_model}} of class \code{SEIR}
##' @include check_arguments.R
##' @export
##' @examples
##' ## Create an SEIR model object.
##' model <- SEIR(
##'   u0 = data.frame(S = 99, E = 0, I = 1, R = 0),
##'   tspan = 1:100,
##'   beta = 0.16,
##'   epsilon = 0.25,
##'   gamma = 0.077
##' )
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
    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments_SEIR())

    ## Check for non-numeric parameters
    check_ldata_arg(nrow(u0), beta, epsilon, gamma)
    beta <- rep(beta, length.out = nrow(u0))
    epsilon <- rep(epsilon, length.out = nrow(u0))
    gamma <- rep(gamma, length.out = nrow(u0))

    ## Arguments seem ok...go on

    G <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1), nrow = 3, ncol = 3,
                dimnames = list(c("S -> beta*S*I/(S+E+I+R) -> E",
                                  "E -> epsilon*E -> I",
                                  "I -> gamma*I -> R"),
                                c("1", "2", "3")))

    S <- matrix(c(-1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1), nrow = 4, ncol = 3,
                dimnames = list(compartments_SEIR(), c("1", "2", "3")))

    ldata <- matrix(as.numeric(c(beta, epsilon, gamma)),
                    nrow  = 3, byrow = TRUE,
                    dimnames = list(c("beta", "epsilon", "gamma")))

    model <- SimInf_model(G      = G,
                          S      = S,
                          E      = select_matrix_SEIR(),
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          u0     = u0)

    methods::as(model, "SEIR")
}

##' Example event data for the \acronym{SEIR} model with cattle herds
##'
##' Dataset containing 466,692 scheduled events for a population of
##' 1,600 cattle herds over 1,460 days (4 years). Demonstrates how
##' demographic and movement events affect SEIR dynamics in a cattle
##' disease context.
##'
##' @details
##' The event data contains three types of scheduled events that affect
##' cattle herds (nodes):
##'
##' \describe{
##'   \item{Exit}{Deaths or removal of cattle from a herd (n =
##'     182,535).  These events decrease the population and affect all
##'     disease compartments proportionally.}
##'   \item{Enter}{Births or introduction of cattle to a herd (n =
##'     182,685).  These events add susceptible cattle to herds,
##'     increasing overall herd size.}
##'   \item{External transfer}{Movement of cattle between herds (n =
##'     101,472).  These events transfer cattle from one herd to
##'     another, potentially facilitating between-herd disease
##'     transmission.}
##' }
##'
##' The \code{select} column in the returned data frame is mapped to
##' the columns of the internal select matrix:
##' \itemize{
##'   \item \code{select = 1} corresponds to \strong{Enter} events,
##'     targeting the Susceptible (S) compartment.
##'   \item \code{select = 2} corresponds to \strong{Exit} and
##'     \strong{External Transfer} events, targeting all compartments.
##' }
##'
##' Events are distributed across all 1,600 herds over the 4-year
##' period, reflecting realistic patterns of cattle demographic change
##' and herd-to-herd movement. The timing and frequency of events can
##' significantly influence disease dynamics simulated by the model.
##'
##' @return A \code{data.frame} with columns:
##'   \describe{
##'     \item{event}{Event type: "exit", "enter", or "extTrans".}
##'     \item{time}{Day when event occurs (1-1460).}
##'     \item{node}{Affected herd identifier (1-1600).}
##'     \item{dest}{Destination herd for external transfer events.}
##'     \item{n}{Number of cattle affected.}
##'     \item{proportion}{0. Not used in this example.}
##'     \item{select}{Model compartment to affect (see
##'       \code{\linkS4class{SimInf_events}}).}
##'     \item{shift}{0. Not used in this example.}
##'   }
##'
##' @seealso
##' \code{\link{u0_SEIR}} for the corresponding initial cattle
##' population, \code{\link{SEIR}} for creating SEIR models with these
##' events, and \code{\linkS4class{SimInf_events}} for event structure
##' details
##'
##' @export
##' @example man/examples/SEIR.R
events_SEIR <- function() {
    utils::data("events_SISe3", package = "SimInf", envir = environment())
    events_SISe3$select[events_SISe3$event == "exit"] <- 2L
    events_SISe3$select[events_SISe3$event == "enter"] <- 1L
    events_SISe3 <- events_SISe3[events_SISe3$event != "intTrans", ]
    events_SISe3$select[events_SISe3$event == "extTrans"] <- 2L
    events_SISe3
}

##' Example Initial population data for the \acronym{SEIR} model
##'
##' Dataset containing the initial number of susceptible, exposed,
##' infected, and recovered cattle across 1,600 herds. Provides
##' realistic population structure for demonstrating SEIR model
##' simulations in a cattle disease epidemiology context.
##'
##' @details
##' This dataset represents initial disease states in a population of
##' 1,600 cattle herds (nodes). Each row represents a single herd
##' (node), derived from the structured cattle population data by
##' adding an exposed compartment to the SIR model structure.
##'
##' The data contains:
##' \describe{
##'   \item{S}{Total susceptible cattle in the herd}
##'   \item{E}{Total exposed cattle (initialized to zero)}
##'   \item{I}{Total infected cattle (initialized to zero)}
##'   \item{R}{Total recovered cattle (initialized to zero)}
##' }
##'
##' The herd size distribution reflects realistic heterogeneity
##' observed in cattle populations, making it suitable for testing
##' disease dynamics with an explicit latent period.
##'
##' @return A \code{data.frame} with 1,600 rows (one per herd) and 4 columns:
##'   \describe{
##'     \item{S}{Number of susceptible cattle in the herd}
##'     \item{E}{Number of exposed cattle in the herd (all zero at start)}
##'     \item{I}{Number of infected cattle in the herd (all zero at start)}
##'     \item{R}{Number of recovered cattle in the herd (all zero at start)}
##'   }
##'
##' @seealso
##' \code{\link{SEIR}} for creating SEIR models with this initial
##' state and \code{\link{events_SEIR}} for associated cattle movement
##' and demographic events
##'
##' @export
##' @example man/examples/SEIR.R
u0_SEIR <- function() {
    u0 <- u0_SIR()
    u0$E <- 0L
    u0[, c("S", "E", "I", "R")]
}
