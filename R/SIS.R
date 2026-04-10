## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
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

##' Class SIS
##'
##' Class to handle the \acronym{SIS} model. This class inherits from
##' \code{\linkS4class{SimInf_model}}, meaning that \acronym{SIS}
##' objects are fully compatible with all generic functions defined
##' for \code{SimInf_model}, such as \code{\link{run}},
##' \code{\link{plot}}, \code{\link{trajectory}}, and
##' \code{\link{prevalence}}.
##'
##' @template SIS-details
##'
##' @seealso
##' \code{\link{SIS}} for creating an \acronym{SIS} model object,
##' \code{\linkS4class{SimInf_model}} for the parent class definition,
##' \code{\link{SIR}} for a model with permanent immunity, and
##' \code{\link{SEIR}} for a model including a latent period.
##' @include SimInf_model.R
##' @export
setClass("SIS", contains = c("SimInf_model"))

##' The compartments in an SIS model
##' @noRd
compartments_SIS <- function() {
    c("S", "I")
}

##' Select matrix for events in the SIS model
##'
##' Internal function returning the 2x2 select matrix (E) that maps
##' SIS compartments (rows) to event types (columns) for event
##' processing.
##'
##' @return A 2x2 numeric matrix with compartments as rows and event
##'     types as columns. Used internally by SimInf_events.
##' @noRd
select_matrix_SIS <- function() {
    matrix(c(1, 0, 1, 1),
           nrow = 2,
           ncol = 2,
           dimnames = list(compartments_SIS(), seq_len(2)))
}

##' Create an SIS model
##'
##' Create an \acronym{SIS} model to be used by the simulation
##' framework.
##'
##' @template SIS-details
##' @details
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of susceptible individuals in each node}
##' \item{I}{The number of infected individuals in each node}
##' }
##'
##' @template u0-param
##' @template tspan-param
##' @template events-param
##' @template beta-param
##' @template gamma-param
##' @return A \code{\link{SimInf_model}} of class \code{SIS}
##' @seealso
##' \code{\linkS4class{SIS}} for the class definition.
##' \code{\link{SIR}}, \code{\link{SEIR}}, \code{\link{SISe}},
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
##' ## Create an SIS model object.
##' model <- SIS(
##'   u0 = data.frame(S = 99, I = 1),
##'   tspan = 1:100,
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Run the SIS model and plot the result.
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
                dimnames = list(c("S -> beta*S*I/(S+I) -> I",
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

##' Example event data for the SIS model with cattle herds
##'
##' Dataset containing 466,692 scheduled events for a population of
##' 1,600 cattle herds over 1,460 days (4 years). Demonstrates how
##' demographic and movement events affect SIS dynamics in a cattle
##' disease context.
##'
##' @details
##' The event data contains three types of scheduled events that
##' affect cattle herds (nodes):
##'
##' \describe{
##'   \item{Exit}{Deaths or removal of cattle from a herd (n =
##'     182,535).  These events decrease the population in both
##'     susceptible and infected compartments.}
##'   \item{Enter}{Births or introduction of cattle to a herd (n =
##'     182,685).  These events add susceptible cattle to herds.}
##'   \item{External transfer}{Movement of cattle between herds (n =
##'     101,472).  These events transfer cattle from one herd to
##'     another, potentially spreading disease across the herd
##'     network. Either susceptible or infected animals may be
##'     transferred.}
##' }
##'
##' The \code{select} column in the returned data frame is mapped to
##' the columns of the internal select matrix:
##' \itemize{
##'   \item \code{select = 1} corresponds to \strong{Enter} events,
##'     targeting the Susceptible (S) compartment.
##'   \item \code{select = 2} corresponds to \strong{Exit} and
##'     \strong{External Transfer} events, targeting all compartments
##'     (S and I).
##' }
##'
##' Events are distributed across all 1,600 herds over the 4-year
##' period, reflecting realistic patterns of cattle demographic change
##' and herd-to-herd movement. In SIS dynamics, these events can
##' introduce disease to previously unaffected herds or remove
##' infected cattle from the system.
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
##' \code{\link{u0_SIS}} for the corresponding initial cattle
##' population, \code{\link{SIS}} for creating SIS models with these
##' events, and \code{\linkS4class{SimInf_events}} for event structure
##' details
##'
##' @export
##' @example man/examples/SIS.R
events_SIS <- function() {
    events_SISe()
}

##' Example initial population data for the SIS model
##'
##' Dataset containing the initial number of susceptible and infected
##' cattle across 1,600 herds. Provides realistic population structure
##' for demonstrating SIS model simulations in a cattle disease
##' epidemiology context.
##'
##' @details
##' This dataset represents initial disease states in a population of
##' 1,600 cattle herds (nodes). Each row represents a single herd
##' (node), derived from the cattle population data by extracting
##' susceptible and infected compartments. The SIS model is
##' appropriate for diseases where recovered individuals do not gain
##' immunity.
##'
##' The data contains:
##' \describe{
##'   \item{S}{Total susceptible cattle in the herd}
##'   \item{I}{Total infected cattle (initialized to zero)}
##' }
##'
##' The herd size distribution reflects realistic heterogeneity
##' observed in cattle populations.
##'
##' @return A \code{data.frame} with 1,600 rows (one per herd) and 2 columns:
##'   \describe{
##'     \item{S}{Number of susceptible cattle in the herd}
##'     \item{I}{Number of infected cattle in the herd (all zero at start)}
##'   }
##'
##' @seealso
##' \code{\link{SIS}} for creating SIS models with this initial state
##' and \code{\link{events_SIS}} for associated cattle movement and
##' demographic events
##'
##' @export
##' @example man/examples/SIS.R
u0_SIS <- function() {
    u0_SISe()
}
