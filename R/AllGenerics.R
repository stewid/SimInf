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
