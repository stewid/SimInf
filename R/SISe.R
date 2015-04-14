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

##' Class \code{"SISe"}
##'
##' Class to handle the SISe \code{\link{siminf_model}}.
##' @name SISe-class
##' @include siminf_model.R
##' @include AllGenerics.R
##' @docType class
##' @keywords classes
##' @export
setClass("SISe", contains = c("siminf_model"))

##' Create a SISe model
##'
##' Create a SISe model to be used by the simulation framework.
##'
##'
##' The argument init must be a \code{data.frame} with the following
##' columns:
##' \describe{
##' \item{id}{Node identifier that uniquely identifies each node. The
##' node identifiers must be zero-based, i.e. the first identifier
##' must be equal to zero.}
##' \item{S}{The number of sucsceptible}
##' \item{I}{The number of infected}
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
##' @param upsilon The response rate from susceptible to infected due
##' to the environmental infectious pressure
##' @param gamma The recover rate from infected to susceptible
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
##' @return \code{SISe}
##' @export
SISe <- function(init,
                 tspan,
                 events  = NULL,
                 phi     = NULL,
                 upsilon = NULL,
                 gamma   = NULL,
                 alpha   = NULL,
                 beta_q1 = NULL,
                 beta_q2 = NULL,
                 beta_q3 = NULL,
                 beta_q4 = NULL,
                 epsilon = NULL)
{
    ## Check init
    if (!all(c("id", "S", "I") %in% names(init))) {
        stop("Missing columns in init")
    }

    init <- init[,c("id", "S", "I")]

    E <- Matrix(c(1, 1,
                  0, 1),
                nrow   = 2,
                ncol   = 2,
                byrow  = TRUE,
                sparse = TRUE)
    E <- as(E, "dgCMatrix")

    S <- new("dgCMatrix")

    G <- Matrix(c(1, 1,
                  1, 1),
                nrow = 2,
                ncol = 2,
                byrow  = TRUE,
                sparse = TRUE)

    N <- Matrix(c(-1,  1,
                   1, -1),
                nrow   = 2,
                ncol   = 2,
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
    if (any(is.null(upsilon),
            is.null(gamma),
            is.null(alpha),
            is.null(beta_q1),
            is.null(beta_q2),
            is.null(beta_q3),
            is.null(beta_q4),
            is.null(epsilon))) {
        stop("Missing parameters to handle the infectious pressure")
    }

    if (!all(is.numeric(upsilon),
             is.numeric(gamma),
             is.numeric(alpha),
             is.numeric(beta_q1),
             is.numeric(beta_q2),
             is.numeric(beta_q3),
             is.numeric(beta_q4),
             is.numeric(epsilon))) {
        stop("Parameters to handle the infectious pressure must be numeric")
    }

    if (!all(identical(length(upsilon), 1L),
             identical(length(gamma), 1L),
             identical(length(alpha), 1L),
             identical(length(epsilon), 1L))) {
        stop("Parameters to handle the infectious pressure must be of length 1")
    }

    upsilon <- rep(upsilon, nrow(init))
    gamma <- rep(gamma, nrow(init))
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
                         upsilon,
                         gamma,
                         alpha,
                         beta_q1,
                         beta_q2,
                         beta_q3,
                         beta_q4,
                         epsilon),
                       nrow  = 9,
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

    return(as(model, "SISe"))
}

##' @rdname run-methods
##' @export
setMethod("run",
          signature(model = "SISe"),
          function(model, threads, seed)
          {
              ## check that siminf_model contains all data structures
              ## required by the siminf solver and that they make sense
              validObject(model);

              .Call(SISe_run, model, threads, seed)
          }
)

##' @rdname susceptible-methods
##' @export
setMethod("susceptible",
          signature("SISe"),
          function(model, ...) {
              as.matrix(model@U[seq(from = 1, to = dim(model@U)[1], by = 2), , drop = FALSE])
          }
)

##' @rdname infected-methods
##' @export
setMethod("infected",
          signature("SISe"),
          function(model, ...) {
              as.matrix(model@U[seq(from = 2, to = dim(model@U)[1], by = 2), , drop = FALSE])
          }
)

##' @rdname prevalence-methods
##' @export
setMethod("prevalence",
          signature("SISe"),
          function(model, ...) {
              I <- colSums(infected(model))
              S <- colSums(susceptible(model))
              I / (S + I)
          }
)

##' plot,-method
##'
##' @param x The \code{model} to plot
##' @param y Unused argument
##' @param ... Additional arguments affecting the plot produced.
##' @name plot-methods
##' @aliases plot plot-methods plot,SISe-method
##' @docType methods
##' @importFrom graphics plot
##' @export
setMethod("plot",
          signature(x = "SISe"),
          function(x, ...)
      {
          savepar <- par(mfrow = c(3,1),
                         mar = c(2,4,1,1),
                         oma = c(2,1,0,0))
          on.exit(par(savepar))

          I <- colSums(infected(x))
          S <- colSums(susceptible(x))

          plot(I / (S + I), type = "l", ylab = "Prevalence")
          plot(I, type = "l", ylab = "Infected")
          plot(S, t = "l", ylab = "Susceptible")

          title(xlab = "Day", outer = TRUE, line = 0)
      }
)
