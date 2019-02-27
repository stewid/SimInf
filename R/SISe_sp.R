## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2019  Stefan Engblom
## Copyright (C) 2015 - 2019  Stefan Widgren
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
##' @template beta-param
##' @param coupling The coupling between neighboring nodes
##' @param distance The distance matrix between neighboring nodes
##' @return \code{SISe_sp}
##' @include check_arguments.R
##' @export
##' @importFrom methods as
##' @importFrom methods is
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
                    distance = NULL)
{
    compartments <- c("S", "I")

    ## Check arguments.

    ## Check u0
    if (!is.data.frame(u0))
        u0 <- as.data.frame(u0)
    if (!all(compartments %in% names(u0)))
        stop("Missing columns in u0")
    u0 <- u0[, compartments, drop = FALSE]

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

    ## Check distance matrix
    if (is.null(distance))
        stop("'distance' is missing")
    if (!is(distance, "dgCMatrix"))
        stop("The 'distance' argument must be of type 'dgCMatrix'")
    if (any(distance < 0))
        stop("All values in the 'distance' matrix must be >= 0")

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
                    nrow = 4, byrow = TRUE)
    ldata <- .Call("SimInf_ldata_sp", ldata, distance, 1L, PACKAGE = "SimInf")

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

    as(model, "SISe_sp")
}
