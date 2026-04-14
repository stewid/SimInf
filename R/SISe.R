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

##' Class SISe
##'
##' Class to handle the \acronym{SISe} \code{\link{SimInf_model}}.
##' @seealso \code{\link{SISe}} for creating an SISe model
##' @include SimInf_model.R
##' @export
setClass("SISe", contains = c("SimInf_model"))

##' The compartments in an SISe model
##' @noRd
compartments_SISe <- function() {
    compartments_SIS()
}

##' Select matrix for events in the SISe model
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

##' Create an SISe model
##'
##' Create an \acronym{SISe} model to be used by the simulation
##' framework.
##'
##' The \acronym{SISe} model contains two compartments: number of
##' susceptible (S) and number of infectious (I)
##' individuals. Additionally, it contains a continuous environmental
##' compartment to model shedding of a pathogen to the
##' environment. Consequently, the model has two state transitions:
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
##' @param gamma A numeric vector with the recovery rate from infected
##'     to susceptible.  Each node can have a different gamma
##'     value. The vector must have length 1 or \code{nrow(u0)}. If
##'     the vector has length 1 but the model contains more nodes, the
##'     value is repeated for all nodes.
##' @template alpha-param
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
