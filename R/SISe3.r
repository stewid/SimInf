## siminf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015  Stefan Engblom
## Copyright (C) 2015  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' Class \code{"SISe3"}
##'
##' Class to handle the SISe3 \code{\link{siminf_model}} model.
##' @name SISe3-class
##' @include siminf_model.r
##' @docType class
##' @keywords classes
##' @export
setClass("SISe3", contains = c("siminf_model"))

##' Coerce a \code{siminf_model} to a \code{SISe3} model
##'
##' @name coerce-siminf_model-SISe3-method
##' @aliases coerce,siminf_model,SISe3-method
##' @docType methods
##' @param from The siminf_model \code{object}
##' @return \code{SISe3} model
##' @keywords methods
setAs(from = "siminf_model", to = "SISe3", def = function(from) {
    return(new("SISe3",
               G            = from@G,
               N            = from@N,
               U            = from@U,
               Nn           = from@Nn,
               data         = from@data,
               sd           = from@sd,
               tspan        = from@tspan,
               u0           = from@u0,
               events       = from@events))
})

##' Create a SISe3 model
##'
##' Create a SISe3 model to be used by the simulation framework.
##'
##'
##' The argument init must be a \code{data.frame} with the following
##' columns:
##' \describe{
##' \item{id}{Node identifier that uniquely identifies each node. The
##' node identifiers must be zero-based, i.e. the first identifier
##' must be equal to zero.}
##' \item{S_age_1}{The number of sucsceptible in age category 1}
##' \item{I_age_1}{The number of infected in age category 1}
##' \item{S_age_2}{The number of sucsceptible in age category 2}
##' \item{I_age_2}{The number of infected in age category 2}
##' \item{S_age_3}{The number of sucsceptible in age category 3}
##' \item{I_age_3}{The number of infected in age category 3}
##' }
##' @param init A \code{data.frame} with the initial state in each
##' node, see details.
##' @param tspan An increasing sequence of points in time where the
##' state of the system is to be returned.
##' @param events a \code{data.frame} with the scheduled events, see
##' \code{\link{siminf_model}}.
##' @param initial_infectious_pressure A numeric vector with the
##' initial environmental infectious pressure in each node. Default
##' NULL which gives 0 in each node.
##' @param response_age_1 The response rate from susceptible to
##' infected for age category 1
##' @param response_age_2 The response rate from susceptible to
##' infected for age category 2
##' @param response_age_3 The response rate from susceptible to
##' infected for age category 3
##' @param recover_age_1 The recover rate from infected to
##' susceptible for age category 1
##' @param recover_age_2 The recover rate from infected to
##' susceptible for age category 2
##' @param recover_age_3 The recover rate from infected to
##' susceptible for age category 3
##' @param alpha The shed rate
##' @param beta_q1 The decay of the environmental infectious pressure
##' in the first quarter of the year.
##' @param beta_q2 The decay of the environmental infectious pressure
##' in the second quarter of the year.
##' @param beta_q3 The decay of the environmental infectious pressure
##' in the third quarter of the year.
##' @param beta_q4 The decay of the environmental infectious pressure
##' in the fourth quarter of the year.
##' @param epsilon The background infectious pressure
##' @return \code{SISe3}
##' @export
SISe3 <- function(init,
                  tspan,
                  events                      = NULL,
                  initial_infectious_pressure = NULL,
                  response_age_1              = NULL,
                  response_age_2              = NULL,
                  response_age_3              = NULL,
                  recover_age_1               = NULL,
                  recover_age_2               = NULL,
                  recover_age_3               = NULL,
                  alpha                       = NULL,
                  beta_q1                     = NULL,
                  beta_q2                     = NULL,
                  beta_q3                     = NULL,
                  beta_q4                     = NULL,
                  epsilon                     = NULL)
{
    ## Check init
    if (!all(c("id",
               "S_age_1",
               "I_age_1",
               "S_age_2",
               "I_age_2",
               "S_age_3",
               "I_age_3") %in% names(init))) {
        stop("Missing columns in init")
    }

    init <- init[,c("id",
                    "S_age_1", "I_age_1",
                    "S_age_2", "I_age_2",
                    "S_age_3", "I_age_3")]

    E <- Matrix(c(1, 0, 0,  1, 0, 0,  2, 0, 0,  1, 0, 0,
                  1, 0, 0,  0, 0, 0,  2, 0, 0,  1, 0, 0,
                  0, 1, 0,  0, 1, 0,  0, 2, 0,  0, 1, 0,
                  0, 1, 0,  0, 0, 0,  0, 2, 0,  0, 1, 0,
                  0, 0, 1,  0, 0, 1,  0, 0, 0,  0, 0, 1,
                  0, 0, 1,  0, 0, 0,  0, 0, 0,  0, 0, 1),
                nrow   = 6,
                ncol   = 12,
                byrow  = TRUE,
                sparse = TRUE)

    G <- Matrix(c(1, 1, 0, 0, 0, 0,
                  1, 1, 0, 0, 0, 0,
                  0, 0, 1, 1, 0, 0,
                  0, 0, 1, 1, 0, 0,
                  0, 0, 0, 0, 1, 1,
                  0, 0, 0, 0, 1, 1),
                nrow   = 6,
                ncol   = 6,
                byrow  = TRUE,
                sparse = TRUE)

    N <- Matrix(c(-1,  1,  0,  0,  0,  0,
                   1, -1,  0,  0,  0,  0,
                   0,  0, -1,  1,  0,  0,
                   0,  0,  1, -1,  0,  0,
                   0,  0,  0,  0, -1,  1,
                   0,  0,  0,  0,  1, -1),
                nrow   = 6,
                ncol   = 6,
                byrow  = TRUE,
                sparse = TRUE)
    N <- as(N, "dgCMatrix")

    ## Check initial infectious pressure
    if (is.null(initial_infectious_pressure))
        initial_infectious_pressure <- rep(0, nrow(init))
    if (!is.numeric(initial_infectious_pressure))
        stop("Invalid 'initial_infectious_pressure': must be numeric vector")
    if (!is.null(dim(initial_infectious_pressure)))
        stop("Invalid 'initial_infectious_pressure': must be numeric vector")
    if (!identical(length(initial_infectious_pressure), nrow(init)))
        stop("Invalid 'initial_infectious_pressure': must be numeric vector with length 'nrow(init)'")
    if (any(initial_infectious_pressure < 0))
        stop("Invalid 'initial_infectious_pressure': must be numeric vector with non-negative values")

    ## Check parameters for response relationship and decay
    if (any(is.null(response_age_1),
            is.null(response_age_2),
            is.null(response_age_3),
            is.null(recover_age_1),
            is.null(recover_age_2),
            is.null(recover_age_3),
            is.null(alpha),
            is.null(beta_q1),
            is.null(beta_q2),
            is.null(beta_q3),
            is.null(beta_q4),
            is.null(epsilon))) {
        stop("Missing parameters to handle the infectious pressure")
    }

    if (!all(is.numeric(response_age_1),
             is.numeric(response_age_2),
             is.numeric(response_age_3),
             is.numeric(recover_age_1),
             is.numeric(recover_age_2),
             is.numeric(recover_age_3),
             is.numeric(alpha),
             is.numeric(beta_q1),
             is.numeric(beta_q2),
             is.numeric(beta_q3),
             is.numeric(beta_q4),
             is.numeric(epsilon))) {
        stop("Parameters to handle the infectious pressure must be numeric")
    }

    if (!all(identical(length(response_age_1), 1L),
             identical(length(response_age_2), 1L),
             identical(length(response_age_3), 1L),
             identical(length(recover_age_1), 1L),
             identical(length(recover_age_2), 1L),
             identical(length(recover_age_3), 1L),
             identical(length(alpha), 1L),
             identical(length(epsilon), 1L))) {
        stop("Parameters to handle the infectious pressure must be of length 1")
    }

    response_age_1 <- rep(response_age_1, nrow(init))
    response_age_2 <- rep(response_age_2, nrow(init))
    response_age_3 <- rep(response_age_3, nrow(init))
    recover_age_1 <- rep(recover_age_1, nrow(init))
    recover_age_2 <- rep(recover_age_2, nrow(init))
    recover_age_3 <- rep(recover_age_3, nrow(init))
    alpha <- rep(alpha, nrow(init))
    epsilon <- rep(epsilon, nrow(init))

    if (identical(length(beta_q1), 1L)) {
        beta_q1 <- rep(beta_q1, nrow(init))
    } else if (!identical(length(beta_q1), nrow(init))) {
        stop("length of beta_q1 must either be 1 or nrow(init)")
    }

    if (identical(length(beta_q2), 1L)) {
        beta_q2 <- rep(beta_q2, nrow(init))
    } else if (!identical(length(beta_q2), nrow(init))) {
        stop("length of beta_q2 must either be 1 or nrow(init)")
    }

    if (identical(length(beta_q3), 1L)) {
        beta_q3 <- rep(beta_q3, nrow(init))
    } else if (!identical(length(beta_q3), nrow(init))) {
        stop("length of beta_q3 must either be 1 or nrow(init)")
    }

    if (identical(length(beta_q4), 1L)) {
        beta_q4 <- rep(beta_q4, nrow(init))
    } else if (!identical(length(beta_q4), nrow(init))) {
        stop("length of beta_q4 must either be 1 or nrow(init)")
    }

    inf_data <- matrix(c(initial_infectious_pressure,
                         response_age_1,
                         response_age_2,
                         response_age_3,
                         recover_age_1,
                         recover_age_2,
                         recover_age_3,
                         alpha,
                         beta_q1,
                         beta_q2,
                         beta_q3,
                         beta_q4,
                         epsilon),
                       nrow = 13,
                       byrow = TRUE)

    storage.mode(inf_data) <- "double"

    model <- siminf_model(G      = G,
                          N      = N,
                          init   = init,
                          E      = E,
                          tspan  = tspan,
                          events = events,
                          data   = inf_data)

    return(as(model, "SISe3"))
}
