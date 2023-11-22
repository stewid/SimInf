## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2022 Stefan Widgren
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

##' A Framework for Data-Driven Stochastic Disease Spread Simulations
##'
##' The SimInf package provides a flexible framework for data-driven
##' spatio-temporal disease spread modeling, designed to efficiently
##' handle population demographics and network data. The framework
##' integrates infection dynamics in each subpopulation as
##' continuous-time Markov chains (CTMC) using the Gillespie
##' stochastic simulation algorithm (SSA) and incorporates available
##' data such as births, deaths or movements as scheduled events. A
##' scheduled event is used to modify the state of a subpopulation at
##' a predefined time-point.
##'
##' The \code{\linkS4class{SimInf_model}} is central and provides the
##' basis for the framework. A \code{\linkS4class{SimInf_model}}
##' object supplies the state-change matrix, the dependency graph, the
##' scheduled events, and the initial state of the system.
##'
##' All predefined models in SimInf have a generating function, with
##' the same name as the model, for example \code{\link{SIR}}.
##'
##' A model can also be created from a model specification using the
##' \code{\link{mparse}} method.
##'
##' After a model is created, a simulation is started with a call to
##' the \code{\link{run}} method and if execution is successful, it
##' returns a modified \code{\linkS4class{SimInf_model}} object with a
##' single stochastic solution trajectory attached to it.
##'
##' SimInf provides several utility functions to inspect simulated
##' data, for example, \code{show}, \code{summary} and \code{plot}.
##' To facilitate custom analysis, it provides the
##' \code{\link{trajectory,SimInf_model-method}} and
##' \code{\link{prevalence}} methods.
##'
##' One of our design goal was to make SimInf extendable and enable
##' usage of the numerical solvers from other R extension packages in
##' order to facilitate complex epidemiological research.  To support
##' this, SimInf has functionality to generate the required C and R
##' code from a model specification, see
##' \code{\link{package_skeleton}}
##' @references
##'
##' \Widgren2019
##' @docType package
##' @aliases SimInf-package
##' @name SimInf
##' @useDynLib SimInf, .registration=TRUE
NULL

##' Unload hook function
##'
##' @param libpath A character string giving the complete path to the
##' package.
##' @noRd
.onUnload <- function(libpath) {
    library.dynam.unload("SimInf", libpath) # nocov
}

##' Example data with spatial distribution of nodes
##'
##' Example data to initialize a population of 1600 nodes and
##' demonstrate various models.
##' @name nodes
##' @docType data
##' @usage data(nodes)
##' @format A \code{data.frame}
##' @keywords dataset
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
##' result <- run(model)
##'
##' ## Determine nodes with one or more infected individuals in the
##' ## trajectory. Extract the 'I' compartment and check for any
##' ## infected individuals in each node.
##' infected <- colSums(trajectory(result, ~ I, format = "matrix")) > 0
##'
##' ## Display infected nodes in 'blue' and non-infected nodes in 'yellow'.
##' data("nodes", package = "SimInf")
##' col <- ifelse(infected, "blue", "yellow")
##' plot(y ~ x, nodes, col = col, pch = 20, cex = 2)
NULL
