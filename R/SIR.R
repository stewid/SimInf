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

##' Class SIR
##'
##' Class to handle the \acronym{SIR} model. This class inherits from
##' \code{\linkS4class{SimInf_model}}, meaning that \acronym{SIR}
##' objects are fully compatible with all generic functions defined
##' for \code{SimInf_model}, such as \code{\link{run}},
##' \code{\link{plot}}, \code{\link{trajectory}}, and
##' \code{\link{prevalence}}.
##'
##' @template SIR-details
##'
##' @seealso
##' \code{\link{SIR}} for creating an \acronym{SIR} model object,
##' \code{\linkS4class{SimInf_model}} for the parent class definition,
##' \code{\link{SEIR}} for a model including a latent period, and
##' \code{\link{SIS}} for a model without immunity.
##' @include SimInf_model.R
##' @export
setClass("SIR", contains = c("SimInf_model"))

##' The compartments in an SIR model
##' @noRd
compartments_SIR <- function() {
    c("S", "I", "R")
}

##' Select matrix for events in the SIR model
##'
##' Internal function returning the 3x4 select matrix (E) that maps
##' SIR compartments (rows) to event types (columns) for event
##' processing.
##'
##' @return A 3x4 numeric matrix with compartments as rows and event
##'     types as columns. Used internally by SimInf_events.
##' @noRd
select_matrix_SIR <- function() {
    matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1),
           nrow = 3,
           ncol = 4,
           dimnames = list(compartments_SIR(), seq_len(4)))
}

##' Create an SIR model
##'
##' Create an \acronym{SIR} model to be used by the simulation
##' framework.
##'
##' @template SIR-details
##' @details
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of susceptible individuals in each node}
##' \item{I}{The number of infected individuals in each node}
##' \item{R}{The number of recovered individuals in each node}
##' }
##'
##' @template u0-param
##' @template tspan-param
##' @template events-param
##' @template beta-param
##' @template gamma-param
##' @return A \code{\link{SimInf_model}} of class \code{SIR}
##' @seealso
##' \code{\linkS4class{SIR}} for the class definition.
##' \code{\link{SEIR}}, \code{\link{SIS}}, \code{\link{SISe}},
##' \code{\link{SISe3}} and \code{\link{SISe_sp}} for other predefined
##' models.  \code{\link{mparse}} for creating custom models.
##' \code{\link{run}} for running the simulation.
##' \code{\link{trajectory}}, \code{\link{prevalence}} and
##' \code{\link{plot,SimInf_model-method}} for post-processing and
##' visualization.
##' @include check_arguments.R
##' @export
##' @examples
##' ## For reproducibility, set the seed.
##' set.seed(22)
##'
##' ## Create an SIR model object.
##' model <- SIR(
##'   u0 = data.frame(S = 99, I = 1, R = 0),
##'   tspan = 1:100,
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Run the SIR model and plot the result.
##' result <- run(model)
##' plot(result)
SIR <- function(u0,
                tspan,
                events = NULL,
                beta   = NULL,
                gamma  = NULL) {
    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments_SIR())

    ## Check for non-numeric parameters
    check_ldata_arg(nrow(u0), beta, gamma)
    beta <- rep(beta, length.out = nrow(u0))
    gamma <- rep(gamma, length.out = nrow(u0))

    ## Arguments seem ok...go on

    G <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2,
                dimnames = list(c("S -> beta*S*I/(S+I+R) -> I",
                                  "I -> gamma*I -> R"),
                                c("1", "2")))

    S <- matrix(c(-1, 1, 0, 0, -1, 1), nrow = 3, ncol = 2,
                dimnames = list(compartments_SIR(), c("1", "2")))

    ldata <- matrix(as.numeric(c(beta, gamma)),
                    nrow  = 2, byrow = TRUE,
                    dimnames = list(c("beta", "gamma")))

    model <- SimInf_model(G      = G,
                          S      = S,
                          E      = select_matrix_SIR(),
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          u0     = u0)

    methods::as(model, "SIR")
}

##' Example event data for the SIR model with cattle herds
##'
##' Dataset containing 466,692 scheduled events for a population of
##' 1,600 cattle herds over 1,460 days (4 years). Demonstrates how
##' demographic and movement events affect SIR dynamics in a cattle
##' disease context.
##'
##' @details
##' The event data contains three types of scheduled events that
##' affect cattle herds (nodes):
##'
##' \describe{
##'   \item{Exit}{Deaths or removal of cattle from a herd (n =
##'     182,535).  These events decrease the population and remove
##'     cattle from the disease system.}
##'   \item{Enter}{Births or introduction of cattle to a herd (n =
##'     182,685).  These events add susceptible cattle to herds,
##'     increasing potential targets for infection.}
##'   \item{External transfer}{Movement of cattle between herds (n =
##'     101,472).  These events transfer cattle from one herd to
##'     another, potentially spreading disease across the herd
##'     network.}
##' }
##'
##' The \code{select} column in the returned data frame is mapped to
##' the columns of the internal select matrix (\code{select_matrix_SIR}):
##' \itemize{
##'   \item \code{select = 1} corresponds to \strong{Enter} events,
##'     targeting the Susceptible (S) compartment.
##'   \item \code{select = 4} corresponds to \strong{Exit} and
##'     \strong{External Transfer} events, targeting all compartments
##'     (S, I, and R).
##' }
##'
##' Events are distributed across all 1,600 herds over the 4-year
##' period. These are synthetic data generated to illustrate how to
##' incorporate scheduled events (such as births, deaths, and
##' movements) into a compartment model in the SimInf framework.
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
##' \code{\link{u0_SIR}} for the corresponding initial cattle
##' population, \code{\link{SIR}} for creating SIR models with these
##' events, and \code{\linkS4class{SimInf_events}} for event structure
##' details
##'
##' @export
##' @example man/examples/SIR.R
events_SIR <- function() {
    utils::data("events_SISe3", package = "SimInf", envir = environment())
    events_SISe3$select[events_SISe3$event == "exit"] <- 4L
    events_SISe3$select[events_SISe3$event == "enter"] <- 1L
    events_SISe3 <- events_SISe3[events_SISe3$event != "intTrans", ]
    events_SISe3$select[events_SISe3$event == "extTrans"] <- 4L
    events_SISe3
}

##' Example initial population data for the SIR model
##'
##' Synthetic dataset containing the initial number of susceptible,
##' infected, and recovered cattle (individuals) across 1,600 cattle
##' herds (nodes).  Provides a heterogeneous population structure for
##' demonstrating SIR model simulations in a compartmental modeling
##' context.
##'
##' @details
##' This dataset represents initial disease states in a synthetic
##' population of 1,600 cattle herds (nodes). Each row represents a
##' single herd (node).
##'
##' The data contains:
##' \describe{
##'   \item{S}{Total susceptible cattle (individuals) in the node}
##'   \item{I}{Total infected cattle (individuals) (initialized to
##'   zero)}
##'   \item{R}{Total recovered cattle (individuals) (initialized to
##'   zero)}
##' }
##'
##' The herd size distribution is synthetically generated to reflect
##' heterogeneity typical of large-scale populations, making it
##' suitable for illustrating how to incorporate scheduled events in
##' the SimInf framework.
##'
##' @return A \code{data.frame} with 1,600 rows (one per node) and 3
##'     columns:
##'     \describe{
##'       \item{S}{Number of susceptible cattle (individuals) in the
##'       herd (node)}
##'       \item{I}{Number of infected cattle (individuals) in the herd
##'       (node) (all zero at start)}
##'       \item{R}{Number of recovered cattle (individuals) in the
##'       herd (node) (all zero at start)}
##'     }
##'
##' @seealso \code{\link{SIR}} for creating SIR models with this
##'     initial state and \code{\link{events_SIR}} for associated
##'     movement and demographic events
##' @export
##' @example man/examples/SIR.R
u0_SIR <- function() {
    utils::data("u0_SISe3", package = "SimInf", envir = environment())
    u0_SISe3$S <- u0_SISe3$S_1 + u0_SISe3$S_2 + u0_SISe3$S_3
    u0_SISe3$I <- u0_SISe3$I_1 + u0_SISe3$I_2 + u0_SISe3$I_3
    u0_SISe3$R <- 0L
    u0_SISe3[, c("S", "I", "R")]
}
