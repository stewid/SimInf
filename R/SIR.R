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

##' Definition of the \acronym{SIR} model
##'
##' Class to handle the \acronym{SIR} \code{\link{SimInf_model}}.
##'
##' The \acronym{SIR} model contains three compartments; number of
##' susceptible (S), number of infectious (I), and number of
##' recovered (R).  Moreover, it has two state transitions,
##' \deqn{S \stackrel{\beta S I / N}{\longrightarrow} I}{
##'   S -- beta S I / N --> I}
##' \deqn{I \stackrel{\gamma I}{\longrightarrow} R}{I -- gamma I --> R}
##' where \eqn{\beta} is the transmission rate, \eqn{\gamma} is the
##' recovery rate, and \eqn{N=S+I+R}.
##' @include SimInf_model.R
##' @export
##' @examples
##' ## Create an SIR model object.
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the SIR model and plot the result.
##' set.seed(22)
##' result <- run(model)
##' plot(result)
setClass("SIR", contains = c("SimInf_model"))

##' The compartments in an SIR model
##' @noRd
compartments_SIR <- function() {
    c("S", "I", "R")
}

##' Select matrix for events in the \acronym{SIR} model
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

##' Create an \acronym{SIR} model
##'
##' Create an \acronym{SIR} model to be used by the simulation
##' framework.
##'
##' The \acronym{SIR} model contains three compartments; number of
##' susceptible (S), number of infectious (I), and number of
##' recovered (R).  Moreover, it has two state transitions,
##' \deqn{S \stackrel{\beta S I / N}{\longrightarrow} I}{
##'   S -- beta S I / N --> I}
##' \deqn{I \stackrel{\gamma I}{\longrightarrow} R}{I -- gamma I --> R}
##' where \eqn{\beta} is the transmission rate, \eqn{\gamma} is the
##' recovery rate, and \eqn{N=S+I+R}.
##'
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of susceptible in each node}
##' \item{I}{The number of infected in each node}
##' \item{R}{The number of recovered in each node}
##' }
##'
##' @template u0-param
##' @template tspan-param
##' @template events-param
##' @template beta-param
##' @template gamma-param
##' @return A \code{\link{SimInf_model}} of class \code{SIR}
##' @include check_arguments.R
##' @export
##' @examples
##' ## Create an SIR model object.
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the SIR model and plot the result.
##' set.seed(22)
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

##' Example Event Data for the \acronym{SIR} Model with Cattle Herds
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
##' Events are distributed across all 1,600 herds over the 4-year
##' period, reflecting realistic patterns of cattle demographic change
##' and herd-to-herd movement in a livestock production system.
##'
##' @return A \code{data.frame} with columns:
##'   \describe{
##'     \item{event}{Event type: "exit", "enter", or "extTrans"}
##'     \item{time}{Day when event occurs (1-1460)}
##'     \item{node}{Affected herd identifier (1-1600)}
##'     \item{dest}{Destination herd for external transfer events}
##'     \item{n}{Number of cattle affected}
##'     \item{select}{Model compartment to affect (see
##'       \code{\linkS4class{SimInf_events}})}
##'   }
##'
##' @seealso
##' \code{\link{u0_SIR}} for the corresponding initial cattle
##' population, \code{\link{SIR}} for creating SIR models with these
##' events, and \code{\linkS4class{SimInf_events}} for event structure
##' details
##'
##' @export
##' @examples
##' ## For reproducibility, call the set.seed() function and specify
##' ## the number of threads to use. To use all available threads,
##' ## remove the set_num_threads() call.
##' set.seed(123)
##' set_num_threads(1)
##'
##' ## Create an 'SIR' model with 1600 cattle herds and initialize
##' ## it to run over 4*365 days. Add one infected animal to the
##' ## first herd to seed the outbreak.
##' u0 <- u0_SIR()
##' u0$I[1] <- 1
##' tspan <- seq(from = 1, to = 4*365, by = 1)
##' model <- SIR(u0     = u0,
##'              tspan  = tspan,
##'              events = events_SIR(),
##'              beta   = 0.16,
##'              gamma  = 0.01)
##'
##' ## Display the number of cattle affected by each event type per day.
##' plot(events(model))
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##'
##' ## Summarize the trajectory.
##' summary(result)
events_SIR <- function() {
    utils::data("events_SISe3", package = "SimInf", envir = environment())
    events_SISe3$select[events_SISe3$event == "exit"] <- 4L
    events_SISe3$select[events_SISe3$event == "enter"] <- 1L
    events_SISe3 <- events_SISe3[events_SISe3$event != "intTrans", ]
    events_SISe3$select[events_SISe3$event == "extTrans"] <- 4L
    events_SISe3
}

##' Example Initial Population Data for the \acronym{SIR} Model
##'
##' Dataset containing the initial number of susceptible, infected,
##' and recovered cattle across 1,600 herds. Provides realistic
##' population structure for demonstrating SIR model simulations in a
##' cattle disease epidemiology context.
##'
##' This dataset represents initial disease states in a population of
##' 1,600 cattle herds. Each node (row) represents a single herd, and
##' the data is derived from the structured \code{u0_SISe3} data by
##' aggregating age-stratified compartments into single S, I, and R
##' compartments for each herd.
##'
##' The aggregated values represent:
##' \describe{
##'   \item{S}{Total susceptible cattle across all age groups in the herd}
##'   \item{I}{Total infected cattle (initialized to zero)}
##'   \item{R}{Total recovered cattle (initialized to zero)}
##' }
##'
##' The herd size distribution reflects realistic heterogeneity
##' observed in cattle populations, making it suitable for testing
##' spatial disease dynamics at the herd level, such as:
##' \itemize{
##'   \item Transmission within and between herds
##'   \item Impact of cattle movement on disease spread
##'   \item Effectiveness of herd-level interventions
##' }
##'
##' @return A \code{data.frame} with 1,600 rows (one per herd) and 3 columns:
##'   \describe{
##'     \item{S}{Number of susceptible cattle in the herd}
##'     \item{I}{Number of infected cattle in the herd (all zero at start)}
##'     \item{R}{Number of recovered cattle in the herd (all zero at start)}
##'   }
##'
##' @seealso
##' \code{\link{SIR}} for creating cattle disease models with this
##' initial state and \code{\link{events_SIR}} for associated cattle
##' movement and demographic events
##'
##' @export
##' @examples
##' \dontrun{
##' ## For reproducibility, call the set.seed() function and specify
##' ## the number of threads to use. To use all available threads,
##' ## remove the set_num_threads() call.
##' set.seed(123)
##' set_num_threads(1)
##'
##' ## Create an 'SIR' model with 1600 cattle herds (nodes) and
##' ## initialize it to run over 4*365 days. Add one infected animal
##' ## to the first herd to seed the outbreak. Define 'tspan' to record
##' ## the state of the system at daily time-points. Load scheduled
##' ## events for the population of nodes with births, deaths and
##' ## between-node movements of individuals.
##' u0 <- u0_SIR()
##' u0$I[1] <- 1
##' model <- SIR(u0     = u0,
##'              tspan  = seq(from = 1, to = 4*365, by = 1),
##'              events = events_SIR(),
##'              beta   = 0.16,
##'              gamma  = 0.01)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##' plot(result)
##'
##' ## Plot the trajectory for the first herd.
##' plot(result, index = 1)
##'
##' ## Summarize trajectory
##' summary(result)
##' }
u0_SIR <- function() {
    utils::data("u0_SISe3", package = "SimInf", envir = environment())
    u0_SISe3$S <- u0_SISe3$S_1 + u0_SISe3$S_2 + u0_SISe3$S_3
    u0_SISe3$I <- u0_SISe3$I_1 + u0_SISe3$I_2 + u0_SISe3$I_3
    u0_SISe3$R <- 0L
    u0_SISe3[, c("S", "I", "R")]
}
