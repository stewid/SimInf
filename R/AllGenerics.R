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

##' Init a \code{SimInf_mparse} object with data
##'
##' A \code{SimInf_mparse} object must be initialised with data to
##' create a \code{SimInf_model} that can be used to simulate from the
##' model.
##' @rdname init-methods
##' @param model The \code{\linkS4class{SimInf_mparse}} object to
##'     initialize.
##' @param u0 A \code{data.frame} (or an object that can be coerced to
##'     a \code{data.frame} with \code{as.data.frame}) with the
##'     initial state in each node.
##' @template tspan-param
##' @param events A \code{data.frame} with the scheduled
##'     events. Default is \code{NULL} i.e. no scheduled events in the
##'     model.
##' @param E Sparse matrix to handle scheduled events, see
##'     \code{\linkS4class{SimInf_events}}. Default is \code{NULL}
##'     i.e. no scheduled events in the model.
##' @param N Sparse matrix to handle scheduled events, see
##'     \code{\linkS4class{SimInf_events}}. Default is \code{NULL}
##'     i.e. no scheduled events in the model.
##' @return a \code{\linkS4class{SimInf_model}} object
##' @template mparse-example
setGeneric("init",
           signature = "model",
           function(model,
                    u0     = NULL,
                    tspan  = NULL,
                    events = NULL,
                    E      = NULL,
                    N      = NULL)
               standardGeneric("init"))

##' Run the SimInf stochastic simulation algorithm
##'
##' @rdname run-methods
##' @param model The siminf model to run.
##' @param threads Number of threads. Default is NULL, i.e. to use all
##'     available processors.
##' @param seed Random number seed. Default is NULL, i.e. the
##'     simulator uses time to seed the random number generator.
##' @param solver Which numerical solver to utilize. Default is Null, i.e.
##'     SSA is the default solver.
##' @return \code{SimInf_model} with result from simulation.
##' @examples
##' ## Create an 'SIR' model with 10 nodes and initialise
##' ## it to run over 100 days.
##' model <- SIR(u0 = data.frame(S = rep(99, 10),
##'                              I = rep(1, 10),
##'                              R = rep(0, 10)),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the model and save the result.
##' result <- run(model, threads = 1, seed = 1)
##'
##' ## Plot the proportion of susceptible, infected and recovered
##' ## individuals.
##' plot(result)
setGeneric("run",
           signature = "model",
           function(model,
                    threads = NULL,
                    seed    = NULL,
                    solver  = NULL)
               standardGeneric("run"))

##' Describe your model in a logical way in R. \code{mparse} creates a
##' \code{\linkS4class{SimInf_mparse}} object with your model
##' definition that is ready to be initialised with data and then
##' \code{\link{run}}.

##' Create a package skeleton for a model depending on SimInf
##'
##' @rdname package_skeleton-methods
##' @param model The \code{model} \code{\linkS4class{SimInf_mparse}}
##'     object with your model to create the package skeleton from.
##' @param name Character string: the package name and directory name
##'     for your package.
##' @param path Path to put the package directory in. Default is '.'
##'     i.e. the current directory.
##' @param author Author of the package.
##' @param email Email of the package maintainer.
##' @param maintainer Maintainer of the package.
##' @param license License of the package. Default is 'GPL-3'.
##' @return invisible \code{NULL}.
##' @export
##' @references Read the \emph{Writing R Extensions} manual for more
##'     details.
##'
##' Once you have created a \emph{source} package you need to install
##' it: see the \emph{R Installation and Administration} manual,
##' \code{\link{INSTALL}} and \code{\link{install.packages}}.
setGeneric("package_skeleton",
           function(model, name = NULL, path = ".", author = NULL,
                    email = NULL, maintainer = NULL,
                    license = "GPL-3") standardGeneric("package_skeleton"))
