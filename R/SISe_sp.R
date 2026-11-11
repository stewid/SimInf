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

##' Class SISe_sp
##'
##' Class to handle the \acronym{SISe_sp} model. This class inherits
##' from \code{\linkS4class{SimInf_model}}, meaning that
##' \acronym{SISe_sp} objects are fully compatible with all generic
##' functions defined for \code{SimInf_model}, such as
##' \code{\link{run}}, \code{\link[=plot,SimInf_model-method]{plot}},
##' \code{\link{trajectory}}, and \code{\link{prevalence}}.
##'
##' @template SISe_sp-details
##'
##' @seealso
##' \code{\link{SISe_sp}} for creating an \acronym{SISe_sp} model
##' object and \code{\linkS4class{SimInf_model}} for the parent class
##' definition.
##' @include SimInf_model.R
##' @export
setClass("SISe_sp", contains = c("SimInf_model"))

##' The compartments in an SISe_sp model
##' @noRd
compartments_SISe_sp <- function() {
    compartments_SIS()
}

##' Select matrix for events in the SISe_sp model
##' @noRd
select_matrix_SISe_sp <- function() {
    select_matrix_SIS()
}

##' Create an SISe_sp model
##'
##' Create a \code{SISe_sp} model to be used by the simulation
##' framework.
##'
##' @template SISe_sp-details
##' @details
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of susceptible}
##' \item{I}{The number of infected}
##' }
##'
##' @template u0-param
##' @template tspan-param
##' @template events-param
##' @template phi-param
##' @param upsilon Indirect transmission rate of the environmental
##'     infectious pressure
##' @param gamma A numeric vector with the recovery rate from infected
##'     to susceptible.  Each node can have a different gamma
##'     value. The vector must have length 1 or \code{nrow(u0)}. If
##'     the vector has length 1 but the model contains more nodes, the
##'     value is repeated for all nodes.
##' @template alpha-param
##' @template beta-end-param
##' @param coupling The coupling between neighboring nodes
##' @param distance The distance matrix between neighboring nodes
##' @return \code{SISe_sp}
##' @seealso \code{\linkS4class{SISe_sp}} for the class definition.
##'     \code{\link{SIR}}, \code{\link{SEIR}}, \code{\link{SIS}},
##'     \code{\link{SISe}} and \code{\link{SISe3_sp}} for other
##'     predefined models.  \code{\link{mparse}} for creating custom
##'     models.  \code{\link{run}} for running the simulation.
##'     \code{\link{trajectory}}, \code{\link{prevalence}} and
##'     \code{\link[=plot,SimInf_model-method]{plot}} for
##'     post-processing and visualization.
##' @include check_arguments.R
##' @export
SISe_sp <- function(u0,
                    tspan,
                    events   = NULL,
                    phi      = NULL,
                    upsilon  = NULL,
                    gamma    = NULL,
                    alpha    = NULL,
                    beta_t1  = NULL,
                    beta_t2  = NULL,
                    beta_t3  = NULL,
                    beta_t4  = NULL,
                    end_t1   = NULL,
                    end_t2   = NULL,
                    end_t3   = NULL,
                    end_t4   = NULL,
                    coupling = NULL,
                    distance = NULL) {
    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_initial_state(u0, compartments_SISe_sp())

    ## Check initial infectious pressure
    if (is.null(phi))
        phi <- 0
    phi <- rep(phi, length.out = nrow(u0))
    check_infectious_pressure_arg(nrow(u0), phi)

    ## Check for non-numeric parameters
    check_gdata_arg(upsilon, gamma, alpha, beta_t1, beta_t2, beta_t3, beta_t4,
                    coupling)

    ## Check interval endpoints
    check_integer_arg(end_t1, end_t2, end_t3, end_t4)
    end_t1 <- rep(end_t1, length.out = nrow(u0))
    end_t2 <- rep(end_t2, length.out = nrow(u0))
    end_t3 <- rep(end_t3, length.out = nrow(u0))
    end_t4 <- rep(end_t4, length.out = nrow(u0))
    check_end_t_arg(nrow(u0), end_t1, end_t2, end_t3, end_t4)

    check_distance_matrix(distance)

    ## Arguments seem ok...go on

    G <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2,
                dimnames = list(c("S -> upsilon*phi*S -> I",
                                  "I -> gamma*I -> S"),
                                c("1", "2")))

    S <- matrix(c(-1,  1, 1, -1), nrow = 2, ncol = 2,
                dimnames = list(compartments_SISe_sp(),
                                c("1", "2")))

    v0 <- matrix(as.numeric(phi), nrow  = 1, byrow = TRUE,
                 dimnames = list("phi"))

    ldata <- matrix(as.numeric(c(end_t1, end_t2, end_t3, end_t4)),
                    nrow = 4, byrow = TRUE,
                    dimnames = list(c("end_t1", "end_t2", "end_t3", "end_t4")))
    ldata <- .Call(SimInf_ldata_sp, ldata, distance, 1L)

    gdata <- as.numeric(c(upsilon, gamma, alpha, beta_t1, beta_t2,
                          beta_t3, beta_t4, coupling))
    names(gdata) <- c("upsilon", "gamma", "alpha", "beta_t1", "beta_t2",
                      "beta_t3", "beta_t4", "coupling")

    model <- SimInf_model(G      = G,
                          S      = S,
                          E      = select_matrix_SISe_sp(),
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0)

    methods::as(model, "SISe_sp")
}
