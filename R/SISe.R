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

##' Definition of the \code{SISe} model
##'
##' Class to handle the SISe \code{\link{SimInf_model}}.
##' @include SimInf_model.R
##' @export
setClass("SISe", contains = c("SimInf_model"))

##' The compartments in an SISe model
##' @noRd
compartments_SISe <- function() {
    compartments_SIS()
}

##' Select matrix for events in the \acronym{SISe} model
##'
##' Internal function returning the 2x2 select matrix (E) that maps
##' SISe compartments (rows) to event types (columns) for event
##' processing.
##'
##' @return A 2x2 numeric matrix with compartments as rows and event
##'     types as columns. Used internally by SimInf_events.
##' @noRd
select_matrix_SISe <- function() {
    select_matrix_SIS()
}

##' Create a SISe model
##'
##' Create an \acronym{SISe} model to be used by the simulation
##' framework.
##'
##' The \acronym{SISe} model contains two compartments; number of
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
##' \item{S}{The number of susceptible in each node}
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
    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments_SISe())

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

    G <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2,
                dimnames = list(c("S -> upsilon*phi*S -> I",
                                  "I -> gamma*I -> S"),
                                c("1", "2")))

    S <- matrix(c(-1,  1, 1, -1), nrow = 2, ncol = 2,
                dimnames = list(compartments_SISe(), c("1", "2")))

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
                          E      = select_matrix_SISe(),
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0)

    methods::as(model, "SISe")
}

##' Example event data for the \acronym{SISe} model with cattle herds
##'
##' Dataset containing 466,692 scheduled events for a population of
##' 1,600 cattle herds over 1,460 days (4 years). Demonstrates how
##' demographic and movement events affect SISe dynamics in a cattle
##' disease context.
##'
##' @details
##' The event data contains three types of scheduled events that affect
##' cattle herds:
##' \describe{
##'   \item{Exit}{Deaths or removal of cattle from a herd (n =
##'     182,535). These events decrease the population in susceptible
##'     and infected compartments.}
##'   \item{Enter}{Births or introduction of cattle to a herd (n =
##'     182,685). These events add susceptible cattle to herds.}
##'   \item{External transfer}{Movement of cattle between herds (n =
##'     101,472). These events transfer cattle from one herd to
##'     another, potentially introducing infected animals.}
##' }
##'
##' Events are distributed across all 1,600 herds over the 4-year
##' period, reflecting realistic patterns of cattle demographic change
##' and herd-to-herd movement.
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
##' \code{\link{u0_SISe}} for the corresponding initial cattle
##' population, \code{\link{SISe}} for creating SISe models with these
##' events and \code{\linkS4class{SimInf_events}} for event structure
##' details
##'
##' @export
##' @example man/examples/SISe.R
events_SISe <- function() {
    utils::data("events_SISe3", package = "SimInf", envir = environment())
    events_SISe3$select[events_SISe3$event == "exit"] <- 2L
    events_SISe3$select[events_SISe3$event == "enter"] <- 1L
    events_SISe3 <- events_SISe3[events_SISe3$event != "intTrans", ]
    events_SISe3$select[events_SISe3$event == "extTrans"] <- 2L
    events_SISe3
}

##' Example initial population data for the \acronym{SISe} model
##'
##' Dataset containing the initial number of susceptible and infected
##' cattle across 1,600 herds, for the environment-based transmission
##' model. Provides realistic population structure for demonstrating
##' SISe model simulations in a cattle disease epidemiology context.
##'
##' @details
##'
##' This dataset represents initial disease states in a population of
##' 1,600 cattle herds (nodes). Each row represents a single herd
##' (node). The SISe model extends the SIS model with an environmental
##' compartment for pathogen shedding, suitable for diseases
##' transmitted through environmental contamination.
##'
##' The data contains:
##' \describe{
##'   \item{S}{Total susceptible cattle in the herd}
##'   \item{I}{Total infected cattle (initialized to zero)}
##' }
##'
##' The herd size distribution reflects realistic heterogeneity
##' observed in cattle populations, making it suitable for testing
##' environmentally- mediated transmission dynamics where pathogen
##' survival in the environment is important.
##'
##' @return A \code{data.frame} with 1,600 rows (one per herd) and 2 columns:
##'   \describe{
##'     \item{S}{Number of susceptible cattle in the herd}
##'     \item{I}{Number of infected cattle in the herd (all zero at start)}
##'   }
##'
##' @seealso
##' \code{\link{SISe}} for creating SISe models with this initial
##' state and \code{\link{events_SISe}} for associated cattle movement
##' and demographic events
##'
##' @export
##' @example man/examples/SISe.R
u0_SISe <- function() {
    u0 <- u0_SIR()
    u0[, c("S", "I")]
}
