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

##' A Framework for Stochastic Disease Spread Simulations
##'
##' @docType package
##' @name SimInf
##' @useDynLib SimInf, .registration=TRUE
NULL

##' Unload hook function
##'
##' @param libpath A character string giving the complete path to the
##' package.
##' @noRd
.onUnload <- function (libpath)
{
    library.dynam.unload("SimInf", libpath)
}

##' Is OpenMP available
##'
##' @return TRUE if SimInf was built with support for OpenMP, else
##'     FALSE.
##' @noRd
have_openmp <- function()
{
    .Call(SimInf_have_openmp)
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
##' result <- run(model, threads = 1, seed = 22)
##'
##' ## Determine nodes with one or more infected individuals in the
##' ## trajectory. Extract the 'I' compartment and check for any
##' ## infected individuals in each node.
##' infected <- colSums(trajectory(result, ~ I, as.is = TRUE)) > 0
##'
##' ## Display infected nodes in 'blue' and non-infected nodes in 'yellow'.
##' data("nodes", package = "SimInf")
##' col <- ifelse(infected, "blue", "yellow")
##' plot(y ~ x, nodes, col = col, pch = 20, cex = 2)
NULL
