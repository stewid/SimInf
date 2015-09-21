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
##' @param beta_t1 The decay of the environmental infectious pressure
##' in the first interval of the year.
##' @param beta_t2 The decay of the environmental infectious pressure
##' in the second interval of the year.
##' @param beta_t3 The decay of the environmental infectious pressure
##' in the third interval of the year.
##' @param beta_t4 The decay of the environmental infectious pressure
##' in the fourth interval of the year.
##' @param epsilon The background infectious pressure
##' @return \code{SISe3}
##' @include check_arguments.R
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
                  beta_t1   = NULL,
                  beta_t2   = NULL,
                  beta_t3   = NULL,
                  beta_t4   = NULL,
                  end_t1    = NULL,
                  end_t2    = NULL,
                  end_t3    = NULL,
                  end_t4    = NULL,
                  epsilon   = NULL)
{
    ## Check arguments.

    ## Check init
    if (!all(c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3") %in% names(init)))
        stop("Missing columns in init")

    ## Check initial infectious pressure
    if (is.null(phi))
        phi <- rep(0, nrow(init))
    check_infectious_pressure_arg(nrow(init), phi)

    ## Check 'gdata' parameters
    check_gdata_arg(upsilon_1, upsilon_2, upsilon_3, gamma_1, gamma_2, gamma_3,
                    alpha, beta_t1, beta_t2, beta_t3, beta_t4, epsilon)

    ## Check interval endpoints
    check_integer_arg(end_t1, end_t2, end_t3, end_t4)
    if (identical(length(end_t1), 1L))
        end_t1 <- rep(end_t1, nrow(init))
    if (identical(length(end_t2), 1L))
        end_t2 <- rep(end_t2, nrow(init))
    if (identical(length(end_t3), 1L))
        end_t3 <- rep(end_t3, nrow(init))
    if (identical(length(end_t4), 1L))
        end_t4 <- rep(end_t4, nrow(init))
    check_end_t_arg(nrow(init), end_t1, end_t2, end_t3, end_t4)

    ## Arguments seems ok...go on

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

    v0 <- matrix(phi, nrow  = 1, byrow = TRUE)
    storage.mode(v0) <- "double"

    ldata <- matrix(c(end_t1, end_t2, end_t3, end_t4),
                    nrow  = 4,
                    byrow = TRUE)
    storage.mode(ldata) <- "double"

    gdata <- c(upsilon_1, upsilon_2, upsilon_3,
               gamma_1, gamma_2, gamma_3,
               alpha,
               beta_t1, beta_t2, beta_t3, beta_t4,
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
          function(model, age = 1:3, i = NULL, by = 1, ...)
          {
              stopifnot(all(age %in% 1:3))

              result <- NULL
              j <- seq(from = 1, to = dim(model@U)[2], by = by)

              if (1 %in% age) {
                  ii <- seq(from = 1, to = dim(model@U)[1], by = 6)
                  if (!is.null(i))
                      ii <- ii[i]
                  result <- as.matrix(model@U[ii, j, drop = FALSE])
              }

              if (2 %in% age) {
                  ii <- seq(from = 3, to = dim(model@U)[1], by = 6)
                  if (!is.null(i))
                      ii <- ii[i]
                  if (is.null(result)) {
                      result <- as.matrix(model@U[ii, j, drop = FALSE])
                  } else {
                      result <- result + as.matrix(model@U[ii, j, drop = FALSE])
                  }
              }

              if (3 %in% age) {
                  ii <- seq(from = 5, to = dim(model@U)[1], by = 6)
                  if (!is.null(i))
                      ii <- ii[i]
                  if (is.null(result)) {
                      result <- as.matrix(model@U[ii, j, drop = FALSE])
                  } else {
                      result <- result + as.matrix(model@U[ii, j, drop = FALSE])
                  }
              }

              result
          }
)

##' @rdname infected-methods
##' @export
setMethod("infected",
          signature("SISe3"),
          function(model, age = 1:3, i = NULL, by = 1, ...)
          {
              stopifnot(all(age %in% 1:3))

              result <- NULL
              j <- seq(from = 1, to = dim(model@U)[2], by = by)

              if (1 %in% age) {
                  ii <- seq(from = 2, to = dim(model@U)[1], by = 6)
                  if (!is.null(i))
                      ii <- ii[i]
                  result <- as.matrix(model@U[ii, j, drop = FALSE])
              }

              if (2 %in% age) {
                  ii <- seq(from = 4, to = dim(model@U)[1], by = 6)
                  if (!is.null(i))
                      ii <- ii[i]
                  if (is.null(result)) {
                      result <- as.matrix(model@U[ii, j, drop = FALSE])
                  } else {
                      result <- result + as.matrix(model@U[ii, j, drop = FALSE])
                  }
              }

              if (3 %in% age) {
                  ii <- seq(from = 6, to = dim(model@U)[1], by = 6)
                  if (!is.null(i))
                      ii <- ii[i]
                  if (is.null(result)) {
                      result <- as.matrix(model@U[ii, j, drop = FALSE])
                  } else {
                      result <- result + as.matrix(model@U[ii, j, drop = FALSE])
                  }
              }

              result
          }
)

##' @rdname prevalence-methods
##' @export
setMethod("prevalence",
          signature("SISe3"),
          function(model, age = 1:3, wnp = FALSE, i = NULL, by = 1, ...)
          {
              I <- infected(model, age = age, i = i, by = by)
              S <- susceptible(model, age = age, i = i, by = by)

              if (identical(wnp, FALSE)) {
                  I <- colSums(I)
                  S <- colSums(S)
              }

              I / (S + I)
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
