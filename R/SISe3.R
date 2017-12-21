## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015 - 2017  Stefan Engblom
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

##' Definition of the \sQuote{SISe3} model
##'
##' Class to handle the SISe3 \code{\link{SimInf_model}} model.
##' @include SimInf_model.R
##' @export
setClass("SISe3", contains = c("SimInf_model"))

##' Create a \sQuote{SISe3} model
##'
##' Create a \sQuote{SISe3} model to be used by the simulation
##' framework.
##'
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
##' @param u0 A \code{data.frame} with the initial state in each
##' node, see details.
##' @template tspan-param
##' @param events a \code{data.frame} with the scheduled events, see
##' \code{\link{SimInf_model}}.
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
##' @template beta-param
##' @param epsilon The background environmental infectious pressure
##' @return \code{SISe3}
##' @include check_arguments.R
##' @export
SISe3 <- function(u0,
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
    compartments <- c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3")

    ## Check arguments.

    ## Check u0
    if (!is.data.frame(u0))
        stop("'u0' must be a data.frame")
    if (!all(compartments %in% names(u0)))
        stop("Missing columns in u0")
    u0 <- u0[, compartments]

    ## Check initial infectious pressure
    if (is.null(phi))
        phi <- rep(0, nrow(u0))
    check_infectious_pressure_arg(nrow(u0), phi)

    ## Check 'gdata' parameters
    check_gdata_arg(upsilon_1, upsilon_2, upsilon_3, gamma_1, gamma_2, gamma_3,
                    alpha, beta_t1, beta_t2, beta_t3, beta_t4, epsilon)

    ## Check interval endpoints
    check_integer_arg(end_t1, end_t2, end_t3, end_t4)
    if (identical(length(end_t1), 1L))
        end_t1 <- rep(end_t1, nrow(u0))
    if (identical(length(end_t2), 1L))
        end_t2 <- rep(end_t2, nrow(u0))
    if (identical(length(end_t3), 1L))
        end_t3 <- rep(end_t3, nrow(u0))
    if (identical(length(end_t4), 1L))
        end_t4 <- rep(end_t4, nrow(u0))
    check_end_t_arg(nrow(u0), end_t1, end_t2, end_t3, end_t4)

    ## Arguments seems ok...go on

    E <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                  1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1),
                nrow = 6, ncol = 6,
                dimnames = list(compartments, c("1", "2", "3", "4", "5", "6")))

    N <- matrix(c(2, 0,
                  2, 0,
                  0, 2,
                  0, 2,
                  0, 0,
                  0, 0),
                nrow   = 6,
                ncol   = 2,
                byrow  = TRUE)
    colnames(N) <- as.character(1:2)
    rownames(N) <- compartments

    G <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                  0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1),
                nrow   = 6, ncol   = 6,
                dimnames = list(c("S_1 -> I_1", "I_1 -> S_1",
                                  "S_2 -> I_2", "I_2 -> S_2",
                                  "S_3 -> I_3", "I_3 -> S_3"),
                                c("1", "2", "3", "4", "5", "6")))

    S <- matrix(c(-1, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0,
                  0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, -1),
                nrow = 6, ncol = 6,
                dimnames = list(compartments, c("1", "2", "3", "4", "5", "6")))

    v0 <- matrix(as.numeric(phi), nrow  = 1, byrow = TRUE)

    ldata <- matrix(c(end_t1, end_t2, end_t3, end_t4), nrow = 4, byrow = TRUE)
    storage.mode(ldata) <- "double"

    gdata <- c(upsilon_1, upsilon_2, upsilon_3,
               gamma_1, gamma_2, gamma_3,
               alpha,
               beta_t1, beta_t2, beta_t3, beta_t4,
               epsilon)
    storage.mode(gdata) <- "double"
    names(gdata) <- c("upsilon_1", "upsilon_2", "upsilon_3",
                      "gamma_1", "gamma_2", "gamma_3",
                      "alpha",
                      "beta_t1", "beta_t2", "beta_t3", "beta_t4",
                      "epsilon")

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

    as(model, "SISe3")
}

##' @rdname plot
##' @aliases plot,SISe3-method
##' @export
##' @importFrom methods callNextMethod
setMethod("plot",
          signature(x = "SISe3"),
          function(x,
                   legend = expression(S[1], I[1], S[2], I[2], S[3], I[3]),
                   col = rep(c("blue", "red"), 3),
                   lty = rep(1:3, each = 2),
                   lwd = 2,
                   ...)
          {
              callNextMethod(x, legend = legend, col = col,
                             lty = lty, lwd = lwd, ...)
          }
)

##' Example data to initialize events for the \sQuote{SISe3} model
##'
##' Example data to initialize scheduled events for a population of
##' 1600 nodes and demonstrate the \code{\linkS4class{SISe3}} model.
##'
##' Example data to initialize scheduled events (see
##' \code{\linkS4class{SimInf_events}}) for a population of 1600 nodes
##' and demonstrate the \code{\linkS4class{SISe3}} model. The dataset
##' contains 783773 events for 1600 nodes distributed over 4 * 365
##' days. The events are divided into three types: \sQuote{Exit}
##' events remove individuals from the population (n = 182535),
##' \sQuote{Enter} events add individuals to the population (n =
##' 182685), sQuote{Internal transfer} events move individuals between
##' compartmens within one node e.g. ageing (n = 317081), and
##' \sQuote{External transfer} events move individuals between nodes
##' in the population (n = 101472). The vignette contains a detailed
##' description of how scheduled events operate on a model.
##' @name events_SISe3
##' @docType data
##' @usage data(events_SISe3)
##' @format A \code{data.frame}
##' @keywords dataset
##' @examples
##' ## Create an 'SISe3' model with 1600 nodes and initialize
##' ## it to run over 4*365 days. Add one infected individual
##' ## to the first node.
##' data("u0_SISe3", package = "SimInf")
##' data("events_SISe3", package = "SimInf")
##' u0_SISe3$I_1[1] <- 1
##' tspan <- seq(from = 1, to = 4*365, by = 1)
##' model <- SISe3(u0 = u0_SISe3, tspan = tspan, events = events_SISe3,
##'                phi = rep(0, nrow(u0_SISe3)), upsilon_1 = 1.8e-2,
##'                upsilon_2 = 1.8e-2, upsilon_3 = 1.8e-2,
##'                gamma_1 = 0.1, gamma_2 = 0.1, gamma_3 = 0.1,
##'                alpha = 1, beta_t1 = 1.0e-1, beta_t2 = 1.0e-1,
##'                beta_t3 = 1.25e-1, beta_t4 = 1.25e-1, end_t1 = 91,
##'                end_t2 = 182, end_t3 = 273, end_t4 = 365, epsilon = 0)
##'
##' ## Display the number of individuals affected by each event type
##' ## per day.
##' plot(events(model))
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model, threads = 1, seed = 3)
##'
##' ## Summarize the trajectory. The summary includes the number of
##' ## events by event type.
##' summary(result)
NULL

##' Example data to initialize the \sQuote{SISe3} model
##'
##' Example data to initialize a population of 1600 nodes and
##' demonstrate the \code{\linkS4class{SISe3}} model.
##'
##' A \code{data.frame} with the number of individuals in the
##' \sQuote{S_1}, \sQuote{S_2}, \sQuote{S_3}, \sQuote{I_1},
##' \sQuote{I_2} and \sQuote{I_3} compartments in 1600 nodes. Note
##' that the \sQuote{I_1}, \sQuote{I_2} and \sQuote{I_3} compartments
##' are zero.
##' @name u0_SISe3
##' @docType data
##' @usage data(u0_SISe3)
##' @format A \code{data.frame}
##' @keywords dataset
##' @examples
##' ## Create an 'SISe3' model with 1600 nodes and initialize it to
##' ## run over 4*365 days and record data at weekly time-points.
##'
##' ## Load the initial population and add ten infected individuals to
##' ## I_1 in the first node.
##' u0 <- u0_SISe3
##' u0$I_1[1] <- 10
##'
##' ## Define 'tspan' to run the simulation over 4*365 and record the
##' ## state of the system at weekly time-points.
##' tspan <- seq(from = 1, to = 4*365, by = 7)
##'
##' ## Load scheduled events for the population of nodes with births,
##' ## deaths and between-node movements of individuals.
##' events <- events_SISe3
##'
##' ## Create a 'SISe3' model
##' model <- SISe3(u0 = u0, tspan = tspan, events = events,
##'                phi = rep(0, nrow(u0)), upsilon_1 = 1.8e-2,
##'                upsilon_2 = 1.8e-2, upsilon_3 = 1.8e-2,
##'                gamma_1 = 0.1, gamma_2 = 0.1, gamma_3 = 0.1,
##'                alpha = 1, beta_t1 = 1.0e-1, beta_t2 = 1.0e-1,
##'                beta_t3 = 1.25e-1, beta_t4 = 1.25e-1, end_t1 = 91,
##'                end_t2 = 182, end_t3 = 273, end_t4 = 365, epsilon = 0)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model, threads = 1, seed = 22)
##'
##' ## Summarize trajectory
##' summary(result)
##'
##' ## Plot the proportion of nodes with at least one infected
##' ## individual.
##' plot(prevalence(result, I_1 + I_2 + I_3 ~ ., "nop"), type = "l")
NULL
