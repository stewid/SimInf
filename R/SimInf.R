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

##' Example data with scheduled events for the \code{SISe3} model
##'
##' Synthetic scheduled events data to demonstrate the \code{SISe3}
##' model. The data contains 783773 events for 1600 nodes over 365 * 4
##' days.
##' @name events_SISe3
##' @docType data
##' @usage data(events_SISe3)
##' @format A \code{data.frame}
##' @keywords dataset
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

##' Example data with spatial distribution of nodes
##'
##' Synthetic data with spatial distribution of 1600 nodes to
##' demonstrate various models.
##' @name nodes
##' @docType data
##' @usage data(nodes)
##' @format A \code{data.frame}
##' @keywords dataset
##' @examples
##' \dontrun{
##' data(nodes)
##' plot(y ~ x, nodes)
##' }
NULL
