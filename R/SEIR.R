## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2017  Stefan Widgren
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

##' Class \code{"SEIR"}
##'
##' Class to handle the SEIR \code{\link{SimInf_model}}.
##' @include SimInf_model.R
##' @include AllGenerics.R
##' @export
setClass("SEIR", contains = c("SimInf_model"))

##' Create a SEIR model
##'
##' Create a SEIR model to be used by the simulation framework.
##'
##'
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S}{The number of sucsceptible in each node}
##' \item{I}{The number of infected in each node}
##' \item{R}{The number of recovered in each node}
##' }
##'
##' @param u0 A \code{data.frame} with the initial state in each node,
##'     see details.
##' @template tspan-param
##' @param events a \code{data.frame} with the scheduled events, see
##'     \code{\link{SimInf_model}}.
##' @param beta The transmission rate from susceptible to exposed.
##' @param epsilon The incubation rate from exposed to infected.
##' @param gamma The recovery rate from infected to recovered.
##' @return \code{SEIR}
##' @include check_arguments.R
##' @export
##' @examples
##' ## Create a SEIR model object.
##' model <- SEIR(u0 = data.frame(S = 99, E = 0, I = 1, R = 0),
##'               tspan = 1:100,
##'               beta = 0.16,
##'               epsilon = 0.25,
##'               gamma = 0.077)
##'
##' ## Run the SEIR model and plot the result.
##' result <- run(model, threads = 1, seed = 3)
##' plot(result)
SEIR <- function(u0,
                 tspan,
                 events  = NULL,
                 beta    = NULL,
                 epsilon = NULL,
                 gamma   = NULL)
{
    compartments <- c("S", "E", "I", "R")

    ## Check arguments.

    ## Check u0
    if (!is.data.frame(u0))
        stop("'u0' must be a data.frame")
    if (!all(compartments %in% names(u0)))
        stop("Missing columns in u0")
    u0 <- u0[, compartments]

    ## Check for non-numeric parameters
    check_gdata_arg(beta, epsilon, gamma)

    ## Arguments seems ok...go on

    E <- Matrix::Matrix(c(1, 1,
                          0, 1,
                          0, 1,
                          0, 1),
                        nrow   = 4,
                        ncol   = 2,
                        byrow  = TRUE,
                        sparse = TRUE)
    E <- methods::as(E, "dgCMatrix")
    colnames(E) <- as.character(1:2)
    rownames(E) <- compartments

    N <- matrix(integer(0), nrow = 0, ncol = 0)

    G <- Matrix::Matrix(c(1, 1, 1,
                          1, 1, 1,
                          1, 1, 1),
                        nrow = 3,
                        ncol = 3,
                        byrow  = TRUE,
                        sparse = TRUE)
    G <- methods::as(G, "dgCMatrix")
    colnames(G) <- as.character(1:3)
    rownames(G) <- c("S -> E", "E -> I", "I -> R")

    S <- Matrix::Matrix(c(-1,  0,  0,
                           1, -1,  0,
                           0,  1, -1,
                           0,  0,  1),
                        nrow   = 4,
                        ncol   = 3,
                        byrow  = TRUE,
                        sparse = TRUE)
    S <- methods::as(S, "dgCMatrix")
    colnames(S) <- as.character(1:3)
    rownames(S) <- compartments

    v0 <- matrix(numeric(0), nrow  = 0, ncol = nrow(u0))
    storage.mode(v0) <- "double"

    ldata <- matrix(numeric(0), nrow = 0, ncol = nrow(u0))
    storage.mode(ldata) <- "double"

    gdata <- c(beta, epsilon, gamma)
    storage.mode(gdata) <- "double"
    names(gdata) <- c("beta", "epsilon", "gamma")

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

    methods::as(model, "SEIR")
}

##' @rdname plot
##' @aliases plot,SEIR-method
##' @export
setMethod("plot",
          signature(x = "SEIR"),
          function(x,
                   col = c("blue", "orange", "red", "darkgreen"),
                   lty = rep(1, 4),
                   lwd = 2,
                   ...)
          {
              methods::callNextMethod(x, col = col, lty = lty, lwd = lwd, ...)
          }
)

##' Example data with scheduled events for the \code{SEIR} model
##'
##' Synthetic scheduled events data to demonstrate the \code{SEIR}
##' model. The data contains 466692 events for 1600 nodes over 365 * 4
##' days.
##' @return A \code{data.frame}
##' @export
events_SEIR <- function() {
    utils::data("events_SISe3", package = "SimInf", envir = environment())
    events_SISe3$select[events_SISe3$event == 0] <- 2
    events_SISe3$select[events_SISe3$event == 1] <- 1
    events_SISe3 <- events_SISe3[events_SISe3$event != 2,]
    events_SISe3$select[events_SISe3$event == 3] <- 2
    events_SISe3
}

##' Example data to initialize the \code{SEIR} model
##'
##' Synthetic init data for 1600 nodes to demonstrate the \code{SEIR}
##' model.
##' @return A \code{data.frame}
##' @export
u0_SEIR <- function() {
    utils::data("u0_SISe3", package = "SimInf", envir = environment())
    u0_SISe3$S <- u0_SISe3$S_1 + u0_SISe3$S_2 + u0_SISe3$S_3
    u0_SISe3$E <- 0
    u0_SISe3$I <- u0_SISe3$I_1 + u0_SISe3$I_2 + u0_SISe3$I_3
    u0_SISe3$R <- 0
    u0_SISe3[, c("S", "E", "I", "R")]
}
