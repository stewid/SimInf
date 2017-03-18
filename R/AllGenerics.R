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
##' @param threads Number of threads. Default is NULL, i.e. to use the
##'     number of available processors.
##' @param seed Random number seed. Default is NULL, i.e. to use a
##'     time-seed.
##' @param U Write the number of individuals in each compartment at
##'     \code{tspan} to the non-zero elements in \code{U}, where
##'     \code{U} is a sparse matrix, \code{dgCMatrix}, with dimension
##'     \eqn{N_n N_c \times} \code{length(tspan)}. Default is
##'     \code{NULL} i.e. to write the number of inidividuals in each
##'     compartment in every node to a dense matrix.
##' @param V Write the real-valued continuous state at \code{tspan} to
##'     the non-zero elements in \code{V}, where \code{V} is a sparse
##'     matrix, \code{dgCMatrix}, with dimension
##'     \eqn{N_n}\code{dim(ldata)[1]} \eqn{\times}
##'     \code{length(tspan)}. Default is \code{NULL} i.e. to write the
##'     real-valued continuous state to a dense matrix.
##' @return \code{SimInf_model} with result from simulation.
##' @examples
##' ## Create a 'SISe' demo model with 1 node and
##' ## initialize it to run over 1000 days.
##' model <- demo_model(nodes = 1, days = 1000, model = "SISe")
##' run(model)
setGeneric("run",
           signature = "model",
           function(model,
                    threads = NULL,
                    seed    = NULL,
                    U       = NULL,
                    V       = NULL)
               standardGeneric("run"))

##' Extract the number of individuals in each compartment
##'
##' The result matrix with the number of individuals in each
##' compartment in every node after running a trajectory with
##' \code{\link{run}}. \code{U[, j]} contains the number of
##' individuals in each compartment at \code{tspan[j]}. \code{U[1:Nc,
##' j]} contains the number of individuals in node 1 at
##' \code{tspan[j]}. \code{U[(Nc + 1):(2 * Nc), j]} contains the
##' number of individuals in node 2 at \code{tspan[j]} etc, where
##' \code{Nc} is the number of compartments in the model. The
##' dimension of the matrix is \eqn{N_n N_c \times}
##' \code{length(tspan)} where \eqn{N_n} is the number of nodes.
##' @rdname U-methods
##' @param model The \code{model} to extract the result matrix from.
##' @return The number of individuals in each compartment
##' @keywords methods
##' @export
##' @examples
##' ## Create an 'SIR' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Run the model
##' result <- run(model, seed = 123)
##'
##' ## Extract the number of individuals in each compartment at the
##' ## time-points in tspan.
##' U(result)
setGeneric("U", function(model) standardGeneric("U"))

##' Set a template for where to write the U result matrix
##'
##' @rdname U_set-methods
##' @param model The \code{model} to set a template for the result
##'     matrix \code{U}.
##' @param value Write the number of individuals in each compartment
##'     at \code{tspan} to the non-zero elements in \code{value},
##'     where \code{value} is a sparse matrix, \code{dgCMatrix}, with
##'     dimension \eqn{N_n N_c \times} \code{length(tspan)}. Default
##'     is \code{NULL} i.e. to write the number of inidividuals in
##'     each compartment in every node to a dense matrix.
##' @keywords methods
##' @export
##' @examples
##' ## Create an 'SIR' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## An example with a sparse U result matrix, which can save a lot
##' ## of memory if the model contains many nodes and time-points, but
##' ## where only a few of the data points are of interest. First
##' ## create a sparse matrix with non-zero entries at the locations
##' ## in U where the number of individuals should be written. Then
##' ## run the model with the sparse matrix as a template for U where
##' ## to write data.
##' m <- Matrix::sparseMatrix(1:18, rep(5:10, each = 3))
##' U(model) <- m
##' result <- run(model, seed = 123)
##'
##' ## Extract the number of individuals in each compartment at the
##' ## time-points in tspan.
##' U(result)
setGeneric("U<-", function(model, value) standardGeneric("U<-"))

##' Extract the continuous state variables
##'
##' The result matrix for the real-valued continuous state. \code{V[,
##' j]} contains the real-valued state of the system at
##' \code{tspan[j]}. The dimension of the matrix is
##' \eqn{N_n}\code{dim(ldata)[1]} \eqn{\times} \code{length(tspan)}.
##' @rdname V-methods
##' @param model The \code{model} to extract the result matrix from.
##' @return The continuous state variables
##' @keywords methods
##' @export
##' @examples
##' ## Create an 'SISe' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6)
##' model <- SISe(u0 = u0, tspan = 1:10,
##'               phi = rep(0, 6),
##'               upsilon = 0.017,
##'               gamma   = 0.1,
##'               alpha   = 1,
##'               beta_t1 = 0.19,
##'               beta_t2 = 0.085,
##'               beta_t3 = 0.075,
##'               beta_t4 = 0.185,
##'               end_t1  = 91,
##'               end_t2  = 182,
##'               end_t3  = 273,
##'               end_t4  = 365,
##'               epsilon = 0.000011)
##'
##' ## Run the model
##' result <- run(model, seed = 123)
##'
##' ## Extract the continuous state variables in each node at the
##' ## time-points in tspan. In the 'SISe' model, V represent the
##' ## environmental infectious pressure phi.
##' V(result)
setGeneric("V", function(model) standardGeneric("V"))

##' Set a template for where to write the V result matrix
##'
##' @rdname V_set-methods
##' @param model The \code{model} to set a template for the result
##'     matrix \code{V}.
##' @param value Write the real-valued continuous state at
##'     \code{tspan} to the non-zero elements in \code{value}, where
##'     \code{value} is a sparse matrix, \code{dgCMatrix}, with
##'     dimension \eqn{N_n}\code{dim(ldata)[1]} \eqn{\times}
##'     \code{length(tspan)}. Default is \code{NULL} i.e. to write the
##'     real-valued continuous state to a dense matrix.
##' @keywords methods
##' @export
##' @examples
##' ## Create an 'SISe' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6)
##' model <- SISe(u0 = u0, tspan = 1:10,
##'               phi = rep(0, 6),
##'               upsilon = 0.017,
##'               gamma   = 0.1,
##'               alpha   = 1,
##'               beta_t1 = 0.19,
##'               beta_t2 = 0.085,
##'               beta_t3 = 0.075,
##'               beta_t4 = 0.185,
##'               end_t1  = 91,
##'               end_t2  = 182,
##'               end_t3  = 273,
##'               end_t4  = 365,
##'               epsilon = 0.000011)
##'
##' ## An example with a sparse V result matrix, which can save a lot
##' ## of memory if the model contains many nodes and time-points, but
##' ## where only a few of the data points are of interest. First
##' ## create a sparse matrix with non-zero entries at the locations
##' ## in V where the continuous state variables should be written. Then
##' ## run the model with the sparse matrix as a template for V where
##' ## to write data.
##' m <- Matrix::sparseMatrix(1:6, 5:10)
##' V(model) <- m
##' result <- run(model, seed = 123)
##'
##' ## Extract the continuous state variables at the time-points in tspan.
##' V(result)
setGeneric("V<-", function(model, value) standardGeneric("V<-"))

##' Susceptible
##'
##' Extracts the number of susceptible.
##' @rdname susceptible-methods
##' @param model The \code{model} to extract the susceptible from
##' @param ... Additional arguments affecting the measure
##' @param age For models with age categories, the age category to
##'     extract.
##' @param i Indices specifying the nodes to include when extracting
##'     the number of susceptible. Default is NULL, which includes all
##'     nodes.
##' @param by The number to increment the sequence of time points
##'     starting from 1. Default is 1, which gives the number of
##'     susceptible at every time point.
##' @keywords methods
##' @export
##' @examples
##' ## Create a 'SISe' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the number of susceptible individuals in each
##' ## node after each time step in the simulation
##' susceptible(result)
##'
##' ## Extract the number of susceptible individuals in the
##' ## first node after each time step in the simulation
##' susceptible(result, i = 1)
##'
##' ## Extract the number of susceptible individuals in the
##' ## first and third node after each time step in the simulation
##' susceptible(result, i = c(1, 3))
##'
##' ## Extract the number of susceptible individuals in the first
##' ## and third node after every other time step in the simulation
##' susceptible(result, i = c(1, 3), by = 2)
##'
##' ## Create a 'SISe3' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe3")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the sum all of susceptible individuals in all age
##' ## categories in each node after each time step in the simulation
##' susceptible(result)
##'
##' ## Extract the number of susceptible individuals in the first age
##' ## category in each node after each time step in the simulation
##' susceptible(result, age = 1)
##'
##' ## Extract the sum of susceptible individuals in the first and
##' ## second age category in each node after each time step in
##' ## the simulation
##' susceptible(result, age = c(1, 2))
##'
##' ## Extract the number of susceptible individuals in the first age
##' ## category in the first and third node after each time step in
##' ## the simulation
##' susceptible(result, i = c(1, 3), age = 1)
setGeneric("susceptible",
           function(model, ...) standardGeneric("susceptible"))

##' Infected
##'
##' Extracts the number of infected
##' @rdname infected-methods
##' @param model The \code{model} to extract the infected from
##' @param ... Additional arguments affecting the measure
##' @param age For models with age categories, the age category to
##' extract.
##' @param i Indices specifying the nodes to include when extracting
##' the number of infected. Default is NULL, which includes all nodes.
##' @param by The number to increment the sequence of time points
##' starting from 1. Default is 1, which gives the number of
##' infected at every time point.
##' @keywords methods
##' @export
##' @examples
##' ## Create a 'SISe' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the number of infected individuals in each
##' ## node after each time step in the simulation
##' infected(result)
##'
##' ## Extract the number of infected individuals in the
##' ## first node after each time step in the simulation
##' infected(result, i = 1)
##'
##' ## Extract the number of infected individuals in the
##' ## first and third node after each time step in the simulation
##' infected(result, i = c(1, 3))
##'
##' ## Extract the number of infected individuals in the first
##' ## and third node after every other time step in the simulation
##' infected(result, i = c(1, 3), by = 2)
##'
##' ## Create a 'SISe3' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe3")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the sum all of infected individuals in all age
##' ## categories in each node after each time step in the simulation
##' infected(result)
##'
##' ## Extract the number of infected individuals in the first age
##' ## category in each node after each time step in the simulation
##' infected(result, age = 1)
##'
##' ## Extract the sum of infected individuals in the first and
##' ## second age category in each node after each time step in
##' ## the simulation
##' infected(result, age = c(1, 2))
##'
##' ## Extract the number of infected individuals in the first age
##' ## category in the first and third node after each time step in
##' ## the simulation
##' infected(result, i = c(1, 3), age = 1)
setGeneric("infected",
           function(model, ...) standardGeneric("infected"))

##' Recovered
##'
##' Extracts the number of recovered
##' @rdname recovered-methods
##' @param model The \code{model} to extract the recovered from
##' @param ... Additional arguments affecting the measure
##' @param i Indices specifying the nodes to include when extracting
##' the number of recovered. Default is NULL, which includes all nodes.
##' @param by The number to increment the sequence of time points
##' starting from 1. Default is 1, which gives the number of
##' recovered at every time point.
##' @keywords methods
##' @export
##' @examples
##' ## Create a 'SIR' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SIR")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the number of recovered individuals in each
##' ## node after each time step in the simulation
##' recovered(result)
##'
##' ## Extract the number of recovered individuals in the
##' ## first node after each time step in the simulation
##' recovered(result, i = 1)
##'
##' ## Extract the number of recovered individuals in the
##' ## first and third node after each time step in the simulation
##' recovered(result, i = c(1, 3))
##'
##' ## Extract the number of recovered individuals in the first
##' ## and third node after every other time step in the simulation
##' recovered(result, i = c(1, 3), by = 2)
setGeneric("recovered",
           function(model, ...) standardGeneric("recovered"))

##' Prevalence
##'
##' Calculate the proportion infected individuals
##' @rdname prevalence-methods
##' @param model The \code{model} to calculated the prevalence from
##' @param ... Additional arguments affecting the measure
##' @param i Indices specifying the nodes to include in the
##' calculation of the prevalence. If \code{wnp = TRUE}, then
##' specifying which nodes to extract prevalence for. Default is NULL,
##' which includes all nodes.
##' @param age For models with age categories, the age category to
##' include in the calculation. Default is that all age categories are
##' included.
##' @param wnp Determine within-node prevalence. Default is FALSE.
##' @param by The number to increment the sequence of time points
##' starting from 1. Default is 1, which gives the prevalence at every
##' time point.
##' @keywords methods
##' @export
##' @examples
##' ## Create a 'SISe' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the prevalence of infected nodes after each time
##' ## step in the simulation
##' prevalence(result)
##'
##' ## Extract the prevalence of infected nodes after each time
##' ## step in the simulation when including only the first,
##' ## second and third node in the population at risk.
##' prevalence(result, i = 1:3)
##'
##' ## Extract the prevalence of infected nodes after every other
##' ## time step in the simulation when including only the first,
##' ## second and third node in the population at risk.
##' prevalence(result, i = 1:3, by = 2)
##'
##' ## Extract the within-node prevalence of infected individuals
##' ## in each node after each time step in the simulation
##' prevalence(result, wnp = TRUE)
##'
##' ## Extract the within-node prevalence of infected individuals
##' ## in the first and third node after each time step in the
##' ## simulation
##' prevalence(result, wnp = TRUE, i = c(1, 3))
##'
##' ## Extract the within-node prevalence of infected individuals
##' ## in the first and third node after every other time step in
##' ## the simulation
##' prevalence(result, wnp = TRUE, i = c(1, 3), by = 2)
##'
##' ## Create a 'SISe3' demo model with 5 nodes and initialize
##' ## it to run over 10 days.
##' model <- demo_model(nodes = 5, days = 10, model = "SISe3")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Extract the prevalence of infected nodes after each time
##' ## step in the simulation
##' prevalence(result)
##'
##' ## Extract the within-node prevalence of infected
##' ## individuals in the third age category after each
##' ## time step in the simulation
##' prevalence(result, wnp = TRUE, age = 3)
setGeneric("prevalence",
           function(model, ...) standardGeneric("prevalence"))
