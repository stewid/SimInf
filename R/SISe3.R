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

##' Class SISe3
##'
##' Class to handle the \acronym{SISe3} model. This class inherits
##' from \code{\linkS4class{SimInf_model}}, meaning that
##' \acronym{SISe3} objects are fully compatible with all generic
##' functions defined for \code{SimInf_model}, such as
##' \code{\link{run}}, \code{\link{plot}}, \code{\link{trajectory}},
##' and \code{\link{prevalence}}.
##'
##' @template SISe3-details
##'
##' @seealso
##' \code{\link{SISe3}} for creating an \acronym{SISe3} model object
##' and \code{\linkS4class{SimInf_model}} for the parent class
##' definition.
##' @include SimInf_model.R
##' @export
setClass("SISe3", contains = c("SimInf_model"))

##' The compartments in an SISe3 model
##' @noRd
compartments_SISe3 <- function() {
    c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3")
}

##' Select matrix for events in the SISe3 model
##'
##' Internal function returning the 6x6 select matrix (E) that maps
##' SISe3 compartments (rows) to event types (columns) for event
##' processing.
##'
##' @return A 6x6 numeric matrix with compartments as rows and event
##'     types as columns. Used internally by SimInf_events.
##' @noRd
select_matrix_SISe3 <- function() {
    matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
             1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1),
           nrow = 6,
           ncol = 6,
           dimnames = list(compartments_SISe3(), seq_len(6)))
}

##' Create a \code{SISe3} model
##'
##' Create a \acronym{SISe3} model to be used by the simulation
##' framework.
##'
##' @template SISe3-details
##' @details
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S_1}{The number of susceptible in age category 1}
##' \item{I_1}{The number of infected in age category 1}
##' \item{S_2}{The number of susceptible in age category 2}
##' \item{I_2}{The number of infected in age category 2}
##' \item{S_3}{The number of susceptible in age category 3}
##' \item{I_3}{The number of infected in age category 3}
##' }
##'
##' @template u0-param
##' @template tspan-param
##' @template events-param
##' @template phi-param
##' @param upsilon_1 Indirect transmission rate of the environmental
##' infectious pressure in age category 1
##' @param upsilon_2 Indirect transmission rate of the environmental
##' infectious pressure in age category 2
##' @param upsilon_3 Indirect transmission rate of the environmental
##' infectious pressure in age category 3
##' @param gamma_1 The recovery rate from infected to susceptible for
##' age category 1
##' @param gamma_2 The recovery rate from infected to susceptible for
##' age category 2
##' @param gamma_3 The recovery rate from infected to susceptible for
##' age category 3
##' @template alpha-param
##' @template beta-end-param
##' @param epsilon The background environmental infectious pressure
##' @return \code{SISe3}
##' @include check_arguments.R
##' @export
SISe3 <- function(u0,
                  tspan,
                  events    = NULL,
                  phi       = NULL,
                  upsilon_1 = NULL,
                  upsilon_2 = NULL,
                  upsilon_3 = NULL,
                  gamma_1   = NULL,
                  gamma_2   = NULL,
                  gamma_3   = NULL,
                  alpha     = NULL,
                  beta_t1   = NULL,
                  beta_t2   = NULL,
                  beta_t3   = NULL,
                  beta_t4   = NULL,
                  end_t1    = NULL,
                  end_t2    = NULL,
                  end_t3    = NULL,
                  end_t4    = NULL,
                  epsilon   = NULL) {

    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments_SISe3())

    ## Check initial infectious pressure
    if (is.null(phi))
        phi <- 0
    phi <- rep(phi, length.out = nrow(u0))
    check_infectious_pressure_arg(nrow(u0), phi)

    ## Check 'gdata' parameters
    check_gdata_arg(upsilon_1, upsilon_2, upsilon_3, gamma_1, gamma_2, gamma_3,
                    alpha, beta_t1, beta_t2, beta_t3, beta_t4, epsilon)

    ## Check interval endpoints
    check_integer_arg(end_t1, end_t2, end_t3, end_t4)
    end_t1 <- rep(end_t1, length.out = nrow(u0))
    end_t2 <- rep(end_t2, length.out = nrow(u0))
    end_t3 <- rep(end_t3, length.out = nrow(u0))
    end_t4 <- rep(end_t4, length.out = nrow(u0))
    check_end_t_arg(nrow(u0), end_t1, end_t2, end_t3, end_t4)

    ## Arguments seem ok...go on

    N <- matrix(c(2, 2, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, -1, 0, -1, 0, -1),
                nrow = 6, ncol = 3,
                dimnames = list(compartments_SISe3(), c("1", "2", "3")))

    G <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                  0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1),
                nrow   = 6, ncol   = 6,
                dimnames = list(c("S_1 -> upsilon_1*phi*S_1 -> I_1",
                                  "I_1 -> gamma_1*I_1 -> S_1",
                                  "S_2 -> upsilon_2*phi*S_2 -> I_2",
                                  "I_2 -> gamma_2*I_2 -> S_2",
                                  "S_3 -> upsilon_3*phi*S_3 -> I_3",
                                  "I_3 -> gamma_3*I_3 -> S_3"),
                                c("1", "2", "3", "4", "5", "6")))

    S <- matrix(c(-1, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0,
                  0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, -1),
                nrow = 6, ncol = 6,
                dimnames = list(compartments_SISe3(),
                                c("1", "2", "3", "4", "5", "6")))

    v0 <- matrix(as.numeric(phi), nrow  = 1, byrow = TRUE,
                 dimnames = list("phi"))

    ldata <- matrix(as.numeric(c(end_t1, end_t2, end_t3, end_t4)),
                    nrow = 4, byrow = TRUE,
                    dimnames = list(c("end_t1", "end_t2", "end_t3", "end_t4")))

    gdata <- as.numeric(c(upsilon_1, upsilon_2, upsilon_3,
                          gamma_1, gamma_2, gamma_3, alpha,
                          beta_t1, beta_t2, beta_t3, beta_t4, epsilon))
    names(gdata) <- c("upsilon_1", "upsilon_2", "upsilon_3",
                      "gamma_1", "gamma_2", "gamma_3", "alpha",
                      "beta_t1", "beta_t2", "beta_t3", "beta_t4", "epsilon")

    model <- SimInf_model(G      = G,
                          S      = S,
                          E      = select_matrix_SISe3(),
                          N      = N,
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0)

    methods::as(model, "SISe3")
}

##' Example event data for the SISe3 model with cattle herds
##'
##' Dataset containing 783,773 scheduled events for a population of
##' 1,600 cattle herds stratified by age over 1,460 days (4
##' years). Demonstrates how demographic, movement, and age-transition
##' events affect SISe3 dynamics in a cattle disease context.
##'
##' @details
##' This dataset contains four types of scheduled events that affect
##' cattle herds (nodes) with age structure:
##'
##' \describe{
##'   \item{Exit}{Deaths or removal of cattle from a herd (n =
##'     182,535). These events remove cattle from susceptible or
##'     infected compartments across age categories.}
##'   \item{Enter}{Births or introduction of cattle to a herd (n =
##'     182,685). These events add susceptible cattle, typically to
##'     the youngest age category.}
##'   \item{Internal transfer}{Age transitions or within-herd
##'     movements (n = 317,081). These events move cattle between age
##'     categories within a herd, reflecting maturation and changing
##'     infection risk with age.}
##'   \item{External transfer}{Movement of cattle between herds (n =
##'     101,472).  These events transfer cattle from one herd to
##'     another across age categories, potentially introducing
##'     infected animals.}
##' }
##'
##' The \code{select} column in the returned data frame is mapped to
##' the columns of the internal select matrix as follows:
##' \itemize{
##'   \item \code{select = 1}: Targets \strong{S_1} (Susceptible, age 1).
##'   \item \code{select = 2}: Targets \strong{S_2} (Susceptible, age 2).
##'   \item \code{select = 3}: Targets \strong{S_3} (Susceptible, age 3).
##'   \item \code{select = 4}: Targets \strong{S_1} and \strong{I_1}
##'     (Susceptible and Infected, age 1).
##'   \item \code{select = 5}: Targets \strong{S_2} and \strong{I_2}
##'     (Susceptible and Infected, age 2).
##'   \item \code{select = 6}: Targets \strong{S_3} and \strong{I_3}
##'     (Susceptible and Infected, age 3).
##' }
##'
##' The \code{shift} column is used for \strong{Internal transfer} events
##' to define the destination compartment. It corresponds to the column
##' index in the internal \code{N} matrix that specifies the transition
##' (e.g., moving from age 1 to age 2).
##'
##' Events are distributed across all 1,600 herds over the 4-year
##' period. These are synthetic data generated to illustrate how to
##' incorporate scheduled events (including births, deaths, movements,
##' and age transitions) into an age-structured compartment model in
##' the SimInf framework. The higher event count compared to
##' non-age-structured models reflects the addition of internal
##' transfer events required for age category transitions.
##'
##' The data contains:
##' \describe{
##'   \item{event}{Event type: "exit", "enter", "intTrans", or "extTrans".}
##'   \item{time}{Day when event occurs (1-1460).}
##'   \item{node}{Affected herd identifier (1-1600).}
##'   \item{dest}{Destination herd for external transfer events, else 0.}
##'   \item{n}{Number of cattle affected.}
##'   \item{select}{Model compartment to affect (see
##'     \code{\linkS4class{SimInf_events}}).}
##'   \item{proportion}{0. Not used in this example.}
##'   \item{shift}{Determines how individuals in internal transfer
##'     events are shifted to enter another compartment.}
##' }
##'
##' @seealso
##' \code{\link{u0_SISe3}} for the corresponding initial cattle
##' population with age structure, \code{\link{SISe3}} for creating
##' SISe3 models with these events and
##' \code{\linkS4class{SimInf_events}} for event structure details
##'
##' @name events_SISe3
##' @docType data
##' @usage data(events_SISe3)
##' @format A \code{data.frame}
##' @keywords dataset
##' @example man/examples/SISe3.R
NULL

##' Example initial population data for the SISe3 model
##'
##' Dataset containing the initial number of susceptible and infected
##' cattle across three age categories in 1,600 herds. Provides
##' realistic population structure for demonstrating SISe3 model
##' simulations in a cattle disease epidemiology context with age
##' structure.
##'
##' @details
##' This dataset represents initial disease states in a population of
##' 1,600 cattle herds (nodes) stratified into three age
##' categories. Each row represents a single herd (node). The SISe3
##' model extends the SISe model with age-structured compartments
##' (S_1, I_1, S_2, I_2, S_3, I_3) and an environmental compartment
##' for pathogen shedding. This is appropriate for diseases where
##' transmission rates or recovery rates differ by age group.
##'
##' The data contains:
##' \describe{
##'   \item{S_1}{Total susceptible cattle in age category 1}
##'   \item{I_1}{Total infected cattle in age category 1 (initialized to zero)}
##'   \item{S_2}{Total susceptible cattle in age category 2}
##'   \item{I_2}{Total infected cattle in age category 2 (initialized to zero)}
##'   \item{S_3}{Total susceptible cattle in age category 3}
##'   \item{I_3}{Total infected cattle in age category 3 (initialized to zero)}
##' }
##'
##' The herd size distribution and age structure reflect realistic
##' heterogeneity observed in cattle populations, making it suitable
##' for testing age-dependent disease dynamics with environmental
##' transmission.
##'
##' @seealso
##' \code{\link{SISe3}} for creating SISe3 models with this initial
##' state and \code{\link{events_SISe3}} for associated cattle
##' movement and demographic events
##' @name u0_SISe3
##' @docType data
##' @usage data(u0_SISe3)
##' @format A \code{data.frame}
##' @keywords dataset
##' @example man/examples/SISe3.R
NULL
