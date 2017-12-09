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

##' Definition of the \code{SIR} model
##'
##' Class to handle the SIR \code{\link{SimInf_model}}.
##' @include SimInf_model.R
##' @export
setClass("SIR", contains = c("SimInf_model"))

##' Create a SIR model
##'
##' Create a SIR model to be used by the simulation framework.
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
##' @param beta The transmission rate from susceptible to infected.
##' @param gamma The recovery rate from infected to recovered.
##' @return \code{SIR}
##' @include check_arguments.R
##' @export
##' @examples
##' ## Create a SIR model object.
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the SIR model and plot the result.
##' result <- run(model, threads = 1, seed = 22)
##' plot(result)
SIR <- function(u0,
                tspan,
                events = NULL,
                beta   = NULL,
                gamma  = NULL)
{
    compartments <- c("S", "I", "R")

    ## Check arguments.

    ## Check u0
    if (!is.data.frame(u0))
        stop("'u0' must be a data.frame")
    if (!all(compartments %in% names(u0)))
        stop("Missing columns in u0")
    u0 <- u0[, compartments]

    ## Check for non-numeric parameters
    check_gdata_arg(beta, gamma)

    ## Arguments seems ok...go on

    E <- Matrix::Matrix(c(1, 1,
                          0, 1,
                          0, 1),
                        nrow   = 3,
                        ncol   = 2,
                        byrow  = TRUE,
                        sparse = TRUE)
    E <- methods::as(E, "dgCMatrix")
    colnames(E) <- as.character(1:2)
    rownames(E) <- compartments

    N <- matrix(integer(0), nrow = 0, ncol = 0)

    G <- Matrix::Matrix(c(1, 1,
                          1, 1),
                        nrow = 2,
                        ncol = 2,
                        byrow  = TRUE,
                        sparse = TRUE)
    G <- methods::as(G, "dgCMatrix")
    colnames(G) <- as.character(1:2)
    rownames(G) <- c("S -> I", "I -> R")

    S <- Matrix::Matrix(c(-1,  0,
                           1, -1,
                           0,  1),
                        nrow   = 3,
                        ncol   = 2,
                        byrow  = TRUE,
                        sparse = TRUE)
    S <- methods::as(S, "dgCMatrix")
    colnames(S) <- as.character(1:2)
    rownames(S) <- compartments

    v0 <- matrix(numeric(0), nrow  = 0, ncol = nrow(u0))
    storage.mode(v0) <- "double"

    ldata <- matrix(numeric(0), nrow = 0, ncol = nrow(u0))
    storage.mode(ldata) <- "double"

    gdata <- c(beta, gamma)
    storage.mode(gdata) <- "double"
    names(gdata) <- c("beta", "gamma")

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

    methods::as(model, "SIR")
}

##' @rdname plot
##' @aliases plot,SIR-method
##' @export
setMethod("plot",
          signature(x = "SIR"),
          function(x,
                   col = c("blue", "red", "darkgreen"),
                   lty = rep(1, 3),
                   lwd = 2,
                   ...)
          {
              methods::callNextMethod(x, col = col, lty = lty,
                                      lwd = lwd, ...)
          }
)

##' Example data to initialize events for the \sQuote{SIR} model
##'
##' Example data to initialize scheduled events for a population of
##' 1600 nodes and demonstrate the \code{\linkS4class{SIR}} model.
##'
##' Example data to initialize scheduled events (see
##' \code{\linkS4class{SimInf_events}}) for a population of 1600 nodes
##' and demonstrate the \code{\linkS4class{SIR}} model. The dataset
##' contains 466692 events for 1600 nodes distributed over 4 * 365
##' days. The events are divided into three types: \sQuote{Exit}
##' events remove individuals from the population (n = 182535),
##' \sQuote{Enter} events add individuals to the population (n =
##' 182685), and \sQuote{External transfer} events move individuals
##' between nodes in the population (n = 101472). The vignette
##' contains a detailed description of how scheduled events operate on
##' a model.
##' @return A \code{data.frame}
##' @export
##' @examples
##' ## Create an 'SIR' model with 1600 nodes and initialize
##' ## it to run over 4*365 days. Add one infected individual
##' ## to the first node.
##' u0 <- u0_SIR()
##' u0$I[1] <- 1
##' tspan <- seq(from = 1, to = 4*365, by = 1)
##' model <- SIR(u0     = u0,
##'              tspan  = tspan,
##'              events = events_SIR(),
##'              beta   = 0.16,
##'              gamma  = 0.077)
##'
##' ## Display the number of individuals affected by each event type
##' ## per day.
##' plot(events(model))
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model, threads = 1, seed = 22)
##'
##' ## Summarize the trajectory. The summary includes the number of
##' ## events by event type.
##' summary(result)
events_SIR <- function() {
    utils::data("events_SISe3", package = "SimInf", envir = environment())
    events_SISe3$select[events_SISe3$event == 0] <- 2
    events_SISe3$select[events_SISe3$event == 1] <- 1
    events_SISe3 <- events_SISe3[events_SISe3$event != 2,]
    events_SISe3$select[events_SISe3$event == 3] <- 2
    events_SISe3
}

##' Example data to initialize the \sQuote{SIR} model
##'
##' Example data to initialize a population of 1600 nodes and
##' demonstrate the \code{\linkS4class{SIR}} model.
##'
##' A \code{data.frame} with the number of individuals in the
##' \sQuote{S}, \sQuote{I} and \sQuote{R} compartments in 1600
##' nodes. Note that the \sQuote{I} and \sQuote{R} compartments are
##' zero.
##' @return A \code{data.frame}
##' @export
##' @examples
##' ## Create an 'SIR' model with 1600 nodes and initialize
##' ## it to run over 4*365 days. Add one infected individual
##' ## to the first node.
##' u0 <- u0_SIR()
##' u0$I[1] <- 1
##' tspan <- seq(from = 1, to = 4*365, by = 1)
##' model <- SIR(u0     = u0,
##'              tspan  = tspan,
##'              events = events_SIR(),
##'              beta   = 0.16,
##'              gamma  = 0.077)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model, threads = 1, seed = 22)
##'
##' ## Summarize trajectory
##' summary(result)
u0_SIR <- function() {
    utils::data("u0_SISe3", package = "SimInf", envir = environment())
    u0_SISe3$S <- u0_SISe3$S_1 + u0_SISe3$S_2 + u0_SISe3$S_3
    u0_SISe3$I <- u0_SISe3$I_1 + u0_SISe3$I_2 + u0_SISe3$I_3
    u0_SISe3$R <- 0
    u0_SISe3[, c("S", "I", "R")]
}
