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
##' @include siminf_model.R
##' @include AllGenerics.R
##' @docType class
##' @keywords classes
##' @export
setClass("SISe3", contains = c("siminf_model"))

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
##' \item{S_1}{The number of sucsceptible in age category 1}
##' \item{I_1}{The number of infected in age category 1}
##' \item{S_2}{The number of sucsceptible in age category 2}
##' \item{I_2}{The number of infected in age category 2}
##' \item{S_3}{The number of sucsceptible in age category 3}
##' \item{I_3}{The number of infected in age category 3}
##' }
##' @param init A \code{data.frame} with the initial state in each
##' node, see details.
##' @param tspan An increasing sequence of points in time where the
##' state of the system is to be returned.
##' @param events a \code{data.frame} with the scheduled events, see
##' \code{\link{siminf_model}}.
##' @param phi A numeric vector with the initial environmental
##' infectious pressure in each node. Default NULL which gives 0 in
##' each node.
##' @param upsilon_1 The response rate from susceptible to
##' infected due to the environmental infectious pressure in age
##' category 1
##' @param upsilon_2 The response rate from susceptible to
##' infected due to the environmental infectious pressure in age
##' category 2
##' @param upsilon_3 The response rate from susceptible to
##' infected due to the environmental infectious pressure in age
##' category 3
##' @param gamma_1 The recover rate from infected to
##' susceptible for age category 1
##' @param gamma_2 The recover rate from infected to
##' susceptible for age category 2
##' @param gamma_3 The recover rate from infected to
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
                  events    = NULL,
                  phi       = NULL,
                  upsilon_1 = NULL,
                  upsilon_2 = NULL,
                  upsilon_3 = NULL,
                  gamma_1   = NULL,
                  gamma_2   = NULL,
                  gamma_3   = NULL,
                  alpha     = NULL,
                  beta_q1   = NULL,
                  beta_q2   = NULL,
                  beta_q3   = NULL,
                  beta_q4   = NULL,
                  epsilon   = NULL)
{
    ## Check init
    if (!all(c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3") %in% names(init))) {
        stop("Missing columns in init")
    }

    init <- init[,c("id",
                    "S_1", "I_1",
                    "S_2", "I_2",
                    "S_3", "I_3")]

    E <- Matrix(c(1, 0, 0, 1, 0, 0,
                  0, 0, 0, 1, 0, 0,
                  0, 1, 0, 0, 1, 0,
                  0, 0, 0, 0, 1, 0,
                  0, 0, 1, 0, 0, 1,
                  0, 0, 0, 0, 0, 1),
                nrow   = 6,
                ncol   = 6,
                byrow  = TRUE,
                sparse = TRUE)

    S <- Matrix(c(2, 0,
                  2, 0,
                  0, 2,
                  0, 2,
                  0, 0,
                  0, 0),
                nrow   = 6,
                ncol   = 2,
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
    if (is.null(phi))
        phi <- rep(0, nrow(init))
    if (!is.numeric(phi))
        stop("Invalid 'phi': must be numeric vector")
    if (!is.null(dim(phi)))
        stop("Invalid 'phi': must be numeric vector")
    if (!identical(length(phi), nrow(init)))
        stop("Invalid 'phi': must be numeric vector with length 'nrow(init)'")
    if (any(phi < 0))
        stop("Invalid 'phi': must be numeric vector with non-negative values")

    ## Check parameters for response relationship and decay
    if (any(is.null(upsilon_1),
            is.null(upsilon_2),
            is.null(upsilon_3),
            is.null(gamma_1),
            is.null(gamma_2),
            is.null(gamma_3),
            is.null(alpha),
            is.null(beta_q1),
            is.null(beta_q2),
            is.null(beta_q3),
            is.null(beta_q4),
            is.null(epsilon))) {
        stop("Missing parameters to handle the infectious pressure")
    }

    if (!all(is.numeric(upsilon_1),
             is.numeric(upsilon_2),
             is.numeric(upsilon_3),
             is.numeric(gamma_1),
             is.numeric(gamma_2),
             is.numeric(gamma_3),
             is.numeric(alpha),
             is.numeric(beta_q1),
             is.numeric(beta_q2),
             is.numeric(beta_q3),
             is.numeric(beta_q4),
             is.numeric(epsilon))) {
        stop("Parameters to handle the infectious pressure must be numeric")
    }

    if (!all(identical(length(upsilon_1), 1L),
             identical(length(upsilon_2), 1L),
             identical(length(upsilon_3), 1L),
             identical(length(gamma_1), 1L),
             identical(length(gamma_2), 1L),
             identical(length(gamma_3), 1L),
             identical(length(alpha), 1L),
             identical(length(epsilon), 1L))) {
        stop("Parameters to handle the infectious pressure must be of length 1")
    }

    upsilon_1 <- rep(upsilon_1, nrow(init))
    upsilon_2 <- rep(upsilon_2, nrow(init))
    upsilon_3 <- rep(upsilon_3, nrow(init))
    gamma_1 <- rep(gamma_1, nrow(init))
    gamma_2 <- rep(gamma_2, nrow(init))
    gamma_3 <- rep(gamma_3, nrow(init))
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

    inf_data <- matrix(c(phi,
                         upsilon_1,
                         upsilon_2,
                         upsilon_3,
                         gamma_1,
                         gamma_2,
                         gamma_3,
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
                          S      = S,
                          tspan  = tspan,
                          events = events,
                          data   = inf_data)

    return(as(model, "SISe3"))
}

##' @rdname run-methods
##' @export
setMethod("run",
          signature(model = "SISe3"),
          function(model, threads, seed)
          {
              ## check that siminf_model contains all data structures
              ## required by the siminf solver and that they make sense
              validObject(model);

              .Call(SISe3_run, model, threads, seed)
          }
)

##' @rdname susceptible-methods
##' @export
setMethod("susceptible",
          signature("SISe3"),
          function(model, age = c("age_1", "age_2", "age_3"), ...) {
              age <- match.arg(age)
              from <- switch(age,
                             age_1 = 1,
                             age_2 = 3,
                             age_3 = 5)
              to = dim(model@U)[1]
              as.matrix(model@U[seq(from = from, to = to, by = 6), , drop = FALSE])
          }
)

##' @rdname infected-methods
##' @export
setMethod("infected",
          signature("SISe3"),
          function(model, age = c("age_1", "age_2", "age_3"), ...) {
              age <- match.arg(age)
              from <- switch(age,
                             age_1 = 2,
                             age_2 = 4,
                             age_3 = 6)
              to = dim(model@U)[1]
              as.matrix(model@U[seq(from = from, to = to, by = 6), , drop = FALSE])
          }
)

##' @name plot-methods
##' @aliases plot plot-methods plot,SISe3-method
##' @docType methods
##' @importFrom graphics plot
##' @export
setMethod("plot",
          signature(x = "SISe3"),
          function(x, t0 = NULL, ...)
      {
          callNextMethod(x,
                         t0 = t0,
                         legend = expression(S[1], I[1], S[2], I[2], S[3], I[3]),
                         ...)
      }
)
