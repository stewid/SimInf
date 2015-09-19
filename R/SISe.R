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
##' @param upsilon Indirect transmission rate of the environmental
##' infectious pressure
##' @param gamma The recovery rate from infected to susceptible
##' @param alpha Shed rate from infected individuals
##' @param beta_t1 The decay of the environmental infectious pressure
##' in the first interval of the year.
##' @param beta_t2 The decay of the environmental infectious pressure
##' in the second interval of the year.
##' @param beta_t3 The decay of the environmental infectious pressure
##' in the third interval of the year.
##' @param beta_t4 The decay of the environmental infectious pressure
##' in the fourth interval of the year.
##' @param epsilon The background infectious pressure
##' @return \code{SISe}
##' @include check_arguments.R
##' @export
SISe <- function(init,
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

    ## Check for missing parameters
    check_null_arg(upsilon, gamma, alpha,
                   beta_t1, beta_t2, beta_t3, beta_t4,
                   end_t1, end_t2, end_t3, end_t4,
                   epsilon)

    # Check for non-numeric parameters
    check_numeric_arg(upsilon, gamma, alpha, beta_t1,
                      beta_t2, beta_t3, beta_t4, epsilon)

    # Check for non-integer parameters
    if (!is.numeric(end_t1))
        stop("'end_t1' must be integer")
    if (!all(is_wholenumber(end_t1)))
        stop("'end_t1' must be integer")
    if (!is.numeric(end_t2))
        stop("'end_t2' must be integer")
    if (!all(is_wholenumber(end_t2)))
        stop("'end_t2' must be integer")
    if (!is.numeric(end_t3))
        stop("'end_t3' must be integer")
    if (!all(is_wholenumber(end_t3)))
        stop("'end_t3' must be integer")
    if (!is.numeric(end_t4))
        stop("'end_t4' must be integer")
    if (!all(is_wholenumber(end_t4)))
        stop("'end_t4' must be integer")

    # Check length of parameters
    if (!identical(length(upsilon), 1L))
        stop("'upsilon' must be of length 1")
    if (!identical(length(gamma), 1L))
        stop("'gamma' must be of length 1")
    if (!identical(length(alpha), 1L))
        stop("'alpha' must be of length 1")
    if (!identical(length(beta_t1), 1L))
        stop("'beta_t1' must be of length 1")
    if (!identical(length(beta_t2), 1L))
        stop("'beta_t2' must be of length 1")
    if (!identical(length(beta_t3), 1L))
        stop("'beta_t3' must be of length 1")
    if (!identical(length(beta_t4), 1L))
        stop("'beta_t4' must be of length 1")
    if (identical(length(end_t1), 1L))
        end_t1 <- rep(end_t1, nrow(init))
    if (!identical(length(end_t1), nrow(init)))
        stop("'end_t1' must be of length 1 or 'nrow(init)'")
    if (identical(length(end_t2), 1L))
        end_t2 <- rep(end_t2, nrow(init))
    if (!identical(length(end_t2), nrow(init)))
        stop("'end_t2' must be of length 1 or 'nrow(init)'")
    if (identical(length(end_t3), 1L))
        end_t3 <- rep(end_t3, nrow(init))
    if (!identical(length(end_t3), nrow(init)))
        stop("'end_t3' must be of length 1 or 'nrow(init)'")
    if (identical(length(end_t4), 1L))
        end_t4 <- rep(end_t4, nrow(init))
    if (!identical(length(end_t4), nrow(init)))
        stop("'end_t4' must be of length 1 or 'nrow(init)'")
    if (!identical(length(epsilon), 1L))
        stop("'epsilon' must be of length 1")

    ## Check interval endpoints
    if (!all(0 <= end_t1))
        stop("'end_t1' must be greater than or equal to '0'")
    if (!all(end_t1 < end_t2))
        stop("'end_t1' must be less than 'end_t2'")
    if (!all(end_t2 < end_t3))
        stop("'end_t2' must be less than 'end_t3'")
    if (!all(end_t3 < 364))
        stop("'end_t3' must be less than '364'")
    if (!all(0 <= end_t4))
        stop("'end_t4' must be greater than or equal to '0'")
    if (!all(end_t4 <= 365))
        stop("'end_t4' must be less than or equal to '365'")
    if (!all((end_t4 < end_t1) | (end_t3 < end_t4)))
        stop("'end_t4' must be less than 'end_t1' or greater than 'end_t3'")

    v0 <- matrix(phi, nrow  = 1, byrow = TRUE)
    storage.mode(v0) <- "double"

    ldata <- matrix(c(end_t1, end_t2, end_t3, end_t4),
                    nrow  = 4,
                    byrow = TRUE)
    storage.mode(ldata) <- "double"

    gdata <- c(upsilon,
               gamma,
               alpha,
               beta_t1,
               beta_t2,
               beta_t3,
               beta_t4,
               epsilon)
    storage.mode(gdata) <- "double"

    model <- siminf_model(G      = G,
                          N      = N,
                          init   = init,
                          E      = E,
                          S      = S,
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          v0     = v0)

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
          function(model, i = NULL, by = 1, ...) {
              ii <- seq(from = 1, to = dim(model@U)[1], by = 2)
              if (!is.null(i))
                  ii <- ii[i]
              j <- seq(from = 1, to = dim(model@U)[2], by = by)
              as.matrix(model@U[ii, j, drop = FALSE])
          }
)

##' @rdname infected-methods
##' @export
setMethod("infected",
          signature("SISe"),
          function(model, i = NULL, by = 1, ...) {
              ii <- seq(from = 2, to = dim(model@U)[1], by = 2)
              if (!is.null(i))
                  ii <- ii[i]
              j <- seq(from = 1, to = dim(model@U)[2], by = by)
              as.matrix(model@U[ii, j, drop = FALSE])
          }
)

##' @rdname prevalence-methods
##' @export
setMethod("prevalence",
          signature("SISe"),
          function(model, wnp = FALSE, i = NULL, by = 1, ...) {
              I <- infected(model = model, i = i, by = by)
              S <- susceptible(model = model, i = i, by = by)

              if (identical(wnp, FALSE)) {
                  I <- colSums(I)
                  S <- colSums(S)
              }

              I / (S + I)
          }
)

##' @name plot-methods
##' @aliases plot plot-methods plot,SISe-method
##' @importFrom graphics plot
##' @export
setMethod("plot",
          signature(x = "SISe"),
          function(x, t0 = NULL, ...)
      {
          callNextMethod(x, t0 = t0, legend = c("S", "I"), ...)
      }
)
