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

##' Definition of the \sQuote{SISe3_sp} model
##'
##' Class to handle the SISe3_sp \code{\link{SimInf_model}} model.
##' @include SimInf_model.R
##' @export
setClass("SISe3_sp", contains = c("SimInf_model"))

##' Create an \code{SISe3_sp} model
##'
##' Create an \code{SISe3_sp} model to be used by the simulation
##' framework.
##'
##' The \code{SISe3_sp} model contains two compartments in three age
##' categories; number of susceptible (S_1, S_2, S_3) and number of
##' infectious (I_1, I_2, I_3). Additionally, it contains an
##' environmental compartment to model shedding of a pathogen to the
##' environment. Moreover, it also includes a spatial coupling of the
##' environmental contamination among proximal nodes to capture
##' between-node spread unrelated to moving infected
##' individuals. Consequently, the model has six state transitions,
##'
##' \deqn{S_1 \stackrel{\upsilon_1 \varphi S_1}{\longrightarrow} I_1}{
##' S_1 -- upsilon_1 phi S_1 --> I_1}
##'
##' \deqn{I_1 \stackrel{\gamma_1 I_1}{\longrightarrow} S_1}{
##' I_1 -- gamma_1 I_1 --> S_1}
##'
##' \deqn{S_2 \stackrel{\upsilon_2 \varphi S_2}{\longrightarrow} I_2}{
##' S_2 -- upsilon_2 phi S_2 --> I_2}
##'
##' \deqn{I_2 \stackrel{\gamma_2 I_2}{\longrightarrow} S_2}{
##' I_2 -- gamma_2 I_2 --> S_2}
##'
##' \deqn{S_3 \stackrel{\upsilon_3 \varphi S_3}{\longrightarrow} I_3}{
##' S_3 -- upsilon_3 phi S_3 --> I_3}
##'
##' \deqn{I_3 \stackrel{\gamma_3 I_3}{\longrightarrow} S_3}{
##' I_3 -- gamma_3 I_3 --> S_3}
##'
##' where the transition rate per unit of time from susceptible to
##' infected is proportional to the concentration of the environmental
##' contamination \eqn{\varphi}{phi} in each node. Moreover, the
##' transition rate from infected to susceptible is the recovery rate
##' \eqn{\gamma_1, \gamma_2, \gamma_3}, measured per individual and
##' per unit of time. Finally, the environmental infectious pressure
##' in each node is evolved by,
##'
##' \deqn{\frac{d \varphi_i(t)}{dt} = \frac{\alpha \left(I_{i,1}(t) +
##' I_{i,2}(t) + I_{i,3}(t)\right)}{N_i(t)} +
##' \sum_k{\frac{\varphi_k(t) N_k(t) - \varphi_i(t) N_i(t)}{N_i(t)}
##' \cdot \frac{D}{d_{ik}}} - \beta(t) \varphi_i(t)}{
##' dphi(t)/dt=
##' alpha (I_1+I_2+I_3)/N+
##' D*sum_k(phi_k*N_k-phi_i*N_i)/(d_ik*N_i)-beta*phi_i}
##'
##' where \eqn{\alpha} is the average shedding rate of the pathogen to
##' the environment per infected individual and \eqn{N = S_1 + S_2 +
##' S_3 + I_1 + I_2 + I_3} the size of the node. Next comes the
##' spatial coupling among proximal nodes, where \eqn{D} is the rate
##' of the local spread and \eqn{d_{ik}} the distance between holdings
##' \eqn{i} and \eqn{k}. The seasonal decay and removal of the
##' pathogen is captured by \eqn{\beta(t)}. The environmental
##' infectious pressure \eqn{\varphi(t)}{phi(t)} in each node is
##' evolved each time unit by the Euler forward method. The value of
##' \eqn{\varphi(t)}{phi(t)} is saved at the time-points specified in
##' \code{tspan}.
##'
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S_1}{The number of sucsceptible in age category 1}
##' \item{I_1}{The number of infected in age category 1}
##' \item{S_2}{The number of sucsceptible in age category 2}
##' \item{I_2}{The number of infected in age category 2}
##' \item{S_3}{The number of sucsceptible in age category 3}
##' \item{I_3}{The number of infected in age category 3}
##' }
##'
##' @template beta-section
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
##' @param alpha Shed rate from infected individuals
##' @template beta-end-param
##' @param coupling The coupling between neighboring nodes
##' @param distance The distance matrix between neighboring nodes
##' @return \code{SISe3_sp}
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
    compartments <- c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3")

    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments)

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

    E <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                  1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1),
                nrow = 6, ncol = 6,
                dimnames = list(compartments, c("1", "2", "3", "4", "5", "6")))

    N <- matrix(c(2, 2, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, -1, 0, -1, 0, -1),
                nrow = 6, ncol = 3,
                dimnames = list(compartments, c("1", "2", "3")))

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
                dimnames = list(compartments, c("1", "2", "3", "4", "5", "6")))

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
                          E      = E,
                          N      = N,
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0)

    methods::as(model, "SISe3_sp")
}
