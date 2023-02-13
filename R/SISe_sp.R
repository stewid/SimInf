## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
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

##' Definition of the \code{SISe_sp} model
##'
##' Class to handle the \code{SISe_sp} \code{\link{SimInf_model}}.
##' @include SimInf_model.R
##' @export
setClass("SISe_sp", contains = c("SimInf_model"))

##' Create a \code{SISe_sp} model
##'
##' Create a \code{SISe_sp} model to be used by the simulation
##' framework.
##'
##' The \code{SISe_sp} model contains two compartments; number of
##' susceptible (S) and number of infectious (I). Additionally, it
##' contains an environmental compartment to model shedding of a
##' pathogen to the environment. Moreover, it also includes a spatial
##' coupling of the environmental contamination among proximal nodes
##' to capture between-node spread unrelated to moving infected
##' individuals. Consequently, the model has two state transitions,
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
##' \deqn{\frac{d \varphi_i(t)}{dt} = \frac{\alpha I_{i}(t)}{N_i(t)} +
##' \sum_k{\frac{\varphi_k(t) N_k(t) - \varphi_i(t) N_i(t)}{N_i(t)}
##' \cdot \frac{D}{d_{ik}}} - \beta(t) \varphi_i(t)}{
##' dphi(t)/dt=
##' alpha I / N +
##' D*sum_k(phi_k*N_k-phi_i*N_i)/(d_ik*N_i)-beta*phi_i}
##'
##' where \eqn{\alpha} is the average shedding rate of the pathogen to
##' the environment per infected individual and \eqn{N = S + I} the
##' size of the node. Next comes the spatial coupling among proximal
##' nodes, where \eqn{D} is the rate of the local spread and
##' \eqn{d_{ik}} the distance between holdings \eqn{i} and
##' \eqn{k}. The seasonal decay and removal of the pathogen is
##' captured by \eqn{\beta(t)}. The environmental infectious pressure
##' \eqn{\varphi(t)}{phi(t)} in each node is evolved each time unit by
##' the Euler forward method. The value of \eqn{\varphi(t)}{phi(t)} is
##' saved at the time-points specified in \code{tspan}.
##'
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of sucsceptible}
##' \item{I}{The number of infected}
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
##' @param coupling The coupling between neighboring nodes
##' @param distance The distance matrix between neighboring nodes
##' @return \code{SISe_sp}
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
    compartments <- c("S", "I")

    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments)

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

    E <- matrix(c(1, 0, 1, 1), nrow = 2, ncol = 2,
                dimnames = list(compartments, c("1", "2")))

    G <- matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2,
                dimnames = list(c("S -> upsilon*phi*S -> I",
                                  "I -> gamma*I -> S"),
                                c("1", "2")))

    S <- matrix(c(-1,  1, 1, -1), nrow = 2, ncol = 2,
                dimnames = list(compartments, c("1", "2")))

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
                          E      = E,
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0)

    methods::as(model, "SISe_sp")
}
