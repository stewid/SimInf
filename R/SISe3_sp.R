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

##' Class SISe3_sp
##'
##' Class to handle the \acronym{SISe3_sp} model. This class inherits
##' from \code{\linkS4class{SimInf_model}}, meaning that
##' \acronym{SISe3_sp} objects are fully compatible with all generic
##' functions defined for \code{SimInf_model}, such as
##' \code{\link{run}}, \code{\link{plot,SimInf_model-method}},
##' \code{\link{trajectory}}, and \code{\link{prevalence}}.
##'
##' @template SISe3_sp-details
##'
##' @seealso
##' \code{\link{SISe3_sp}} for creating an \acronym{SISe3_sp} model
##' object and \code{\linkS4class{SimInf_model}} for the parent class
##' definition.
##' @include SimInf_model.R
##' @export
setClass("SISe3_sp", contains = c("SimInf_model"))

##' The compartments in an SISe3_sp model
##' @noRd
compartments_SISe3_sp <- function() {
    compartments_SISe3()
}

##' Select matrix for events in the SISe3_sp model
##'
##' Internal function returning the 6x6 select matrix (E) that maps
##' SISe3_sp compartments (rows) to event types (columns) for event
##' processing.
##'
##' @return A 6x6 numeric matrix with compartments as rows and event
##'     types as columns. Used internally by SimInf_events.
##' @noRd
select_matrix_SISe3_sp <- function() {
    select_matrix_SISe3()
}

##' Create an SISe3_sp model
##'
##' Create an \code{SISe3_sp} model to be used by the simulation
##' framework.
##'
##' @template SISe3_sp-details
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
##' @param coupling The coupling between neighboring nodes
##' @param distance The distance matrix between neighboring nodes
##' @return \code{SISe3_sp}
##' @seealso
##' \code{\linkS4class{SISe3_sp}} for the class definition.
##' \code{\link{SIR}}, \code{\link{SEIR}}, \code{\link{SIS}},
##' \code{\link{SISe3}} and \code{\link{SISe_sp}} for other predefined
##' models.  \code{\link{mparse}} for creating custom models.
##' \code{\link{run}} for running the simulation.
##' \code{\link{trajectory}}, \code{\link{prevalence}} and
##' \code{\link{plot,SimInf_model-method}} for post-processing and
##' visualization.
##' @include check_arguments.R
##' @export
SISe3_sp <- function(u0,
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
                     distance  = NULL,
                     coupling  = NULL) {
    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_initial_state(u0, compartments_SISe3_sp())

    ## Check initial infectious pressure
    if (is.null(phi))
        phi <- 0
    phi <- rep(phi, length.out = nrow(u0))
    check_infectious_pressure_arg(nrow(u0), phi)

    ## Check 'gdata' parameters
    check_gdata_arg(upsilon_1, upsilon_2, upsilon_3, gamma_1, gamma_2, gamma_3,
                    alpha, beta_t1, beta_t2, beta_t3, beta_t4, coupling)

    ## Check interval endpoints
    check_integer_arg(end_t1, end_t2, end_t3, end_t4)
    end_t1 <- rep(end_t1, length.out = nrow(u0))
    end_t2 <- rep(end_t2, length.out = nrow(u0))
    end_t3 <- rep(end_t3, length.out = nrow(u0))
    end_t4 <- rep(end_t4, length.out = nrow(u0))
    check_end_t_arg(nrow(u0), end_t1, end_t2, end_t3, end_t4)

    check_distance_matrix(distance)

    ## Arguments seem ok...go on

    N <- matrix(c(2, 2, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, -1, 0, -1, 0, -1),
                nrow = 6, ncol = 3,
                dimnames = list(compartments_SISe3_sp(),
                                c("1", "2", "3")))

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
                dimnames = list(compartments_SISe3_sp(),
                                c("1", "2", "3", "4", "5", "6")))

    v0 <- matrix(as.numeric(phi), nrow  = 1, byrow = TRUE,
                 dimnames = list("phi"))

    ldata <- matrix(as.numeric(c(end_t1, end_t2, end_t3, end_t4)),
                    nrow = 4, byrow = TRUE,
                    dimnames = list(c("end_t1", "end_t2", "end_t3", "end_t4")))
    ldata <- .Call(SimInf_ldata_sp, ldata, distance, 1L)

    gdata <- as.numeric(c(upsilon_1, upsilon_2, upsilon_3,
                          gamma_1, gamma_2, gamma_3, alpha,
                          beta_t1, beta_t2, beta_t3, beta_t4, coupling))
    names(gdata) <- c("upsilon_1", "upsilon_2", "upsilon_3",
                      "gamma_1", "gamma_2", "gamma_3", "alpha",
                      "beta_t1", "beta_t2", "beta_t3", "beta_t4", "coupling")

    model <- SimInf_model(G      = G,
                          S      = S,
                          E      = select_matrix_SISe3_sp(),
                          N      = N,
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0)

    methods::as(model, "SISe3_sp")
}
