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
##' @keywords internal
.onUnload <- function (libpath)
{
    library.dynam.unload("SimInf", libpath)
}

##' Is OpenMP available
##'
##' @return TRUE if SimInf was built with support for OpenMP, else
##'     FALSE.
##' @keywords internal
have_openmp <- function()
{
    .Call(SimInf_have_openmp)
}

##' Scheduled events example data
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

##' Example data to initialize a model
##'
##' Synthetic init data for 1600 nodes to demonstrate the \code{SISe3}
##' model.
##' @name u0_SISe3
##' @docType data
##' @usage data(u0_SISe3)
##' @format A \code{data.frame}
##' @keywords dataset
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
