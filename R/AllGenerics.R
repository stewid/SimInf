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

##' Extract the C code from an \code{mparse} object
##'
##' @rdname C_code-methods
##' @param model The \code{mparse} object to extract the C code from.
##' @param pkg Character vector. If the C could should be used in a
##'     package named \code{pkg}, the function modifies the C code to
##'     facilitate adding the code to the package. Default is to not
##'     use this argument and return the C code unmodified.
##' @return Character vector with C code for the model.
##' @export
##' @examples
##' ## Use the model parser to create a 'SimInf_mparse' object that
##' ## expresses an SIR model, where 'b' is the transmission rate and
##' ## 'g' is the recovery rate.
##' m <- mparse(c("S -> b*S*I/(S+I+R) -> I", "I -> g*I -> R"),
##'             c("S", "I", "R"), b = 0.16, g = 0.077)
##'
##' ## View the C code.
##' C_code(m)
##'
##' ## Modify the C code for a package named "XYZ"
##' C_code(m, "XYZ")
setGeneric("C_code", function(model, pkg) standardGeneric("C_code"))

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

##' Extract the number of individuals in each compartment
##'
##' The number of individuals in each compartment in every node after
##' running a trajectory with \code{\link{run}}.
##'
##' Description of the layout of the matrix that is returned if
##' \code{as.is = TRUE}. \code{U[, j]} contains the number of
##' individuals in each compartment at \code{tspan[j]}. \code{U[1:Nc,
##' j]} contains the number of individuals in node 1 at
##' \code{tspan[j]}. \code{U[(Nc + 1):(2 * Nc), j]} contains the
##' number of individuals in node 2 at \code{tspan[j]} etc, where
##' \code{Nc} is the number of compartments in the model. The
##' dimension of the matrix is \eqn{N_n N_c \times}
##' \code{length(tspan)} where \eqn{N_n} is the number of nodes.
##' @rdname U-methods
##' @param model the \code{model} to extract the result from.
##' @param as.is the default (\code{as.is = FALSE}) is to generate a
##'     \code{data.frame} with one row per node and time-step with the
##'     number of individuals in each compartment. Using \code{as.is =
##'     TRUE} returns the result as a matrix, which is the internal
##'     format (see \sQuote{Details}).
##' @return The number of individuals in each compartment
##' @export
##' @examples
##' ## Create an 'SIR' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Run the model
##' result <- run(model, threads = 1, seed = 22)
##'
##' ## Extract the number of individuals in each compartment at the
##' ## time-points in tspan.
##' U(result)
setGeneric("U", function(model, as.is = FALSE) standardGeneric("U"))

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
##' result <- run(model, threads = 1, seed = 22)
##'
##' ## Extract the number of individuals in each compartment at the
##' ## time-points in tspan.
##' U(result)
setGeneric("U<-", function(model, value) standardGeneric("U<-"))

##' Extract the continuous state variables
##'
##' The continuous state variables in every node after running a
##' trajectory with \code{\link{run}}.
##'
##' Description of the layout of the matrix that is returned if
##' \code{as.is = TRUE}. The result matrix for the real-valued
##' continuous state. \code{V[, j]} contains the real-valued state of
##' the system at \code{tspan[j]}. The dimension of the matrix is
##' \eqn{N_n}\code{dim(ldata)[1]} \eqn{\times} \code{length(tspan)}.
##' @rdname V-methods
##' @param model The \code{model} to extract the result matrix from.
##' @param as.is the default (\code{as.is = FALSE}) is to generate a
##'     \code{data.frame} with one row per node and time-step with the
##'     value of continuous state variables. Using \code{as.is = TRUE}
##'     returns the result as a matrix, which is the internal format
##'     (see \sQuote{Details}).
##' @return The continuous state variables
##' @export
##' @examples
##' ## Create an 'SISe' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6)
##' model <- SISe(u0 = u0, tspan = 1:10, phi = rep(0, 6),
##'     upsilon = 0.02, gamma = 0.1, alpha = 1, epsilon = 1.1e-5,
##'     beta_t1 = 0.15, beta_t2 = 0.15, beta_t3 = 0.15, beta_t4 = 0.15,
##'     end_t1 = 91, end_t2 = 182, end_t3 = 273, end_t4 = 365)
##'
##' ## Run the model
##' result <- run(model, threads = 1, seed = 7)
##'
##' ## Extract the continuous state variables in each node at the
##' ## time-points in tspan. In the 'SISe' model, V represent the
##' ## environmental infectious pressure phi.
##' V(result)
setGeneric("V", function(model, as.is = FALSE) standardGeneric("V"))

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
##' @export
##' @examples
##' ## Create an 'SISe' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6)
##' model <- SISe(u0 = u0, tspan = 1:10, phi = rep(0, 6),
##'     upsilon = 0.02, gamma = 0.1, alpha = 1, epsilon = 1.1e-5,
##'     beta_t1 = 0.15, beta_t2 = 0.15, beta_t3 = 0.15, beta_t4 = 0.15,
##'     end_t1 = 91, end_t2 = 182, end_t3 = 273, end_t4 = 365)
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
##' result <- run(model, threads = 1, seed = 7)
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
##' @export
##' @examples
##' ## Create an 'SIR' model with 5 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = rep(99, 5), I = rep(1, 5), R = rep(0, 5))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
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
##' @export
##' @examples
##' ## Create an 'SIR' model with 5 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = rep(99, 5), I = rep(1, 5), R = rep(0, 5))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Run the model and save the result
##' result <- run(model, threads = 1, seed = 1)
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
##' ## first and fifth node after each time step in the simulation
##' infected(result, i = c(1, 5))
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
##' @export
##' @examples
##' ## Create an 'SIR' model with 5 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = rep(99, 5), I = rep(1, 5), R = rep(0, 5))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Run the model and save the result
##' result <- run(model, threads = 1, seed = 1)
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
setGeneric("recovered",
           function(model, ...) standardGeneric("recovered"))

##' Prevalence
##'
##' Calculate the proportion of individuals with disease, or the
##' proportion of nodes with individuals with disease, or the
##' proportion of individuals with disease in each node.
##' @rdname prevalence-methods
##' @param model The \code{model} to calculated the prevalence from.
##' @param type The type of prevalence measure to calculate:
##'     \code{'pop'} (default) calcalates the proportion of the
##'     individuals in the population that have disease (model
##'     specific) at each time point in \code{tspan}, \code{'bnp'}
##'     calculates the between-node prevalence, and \code{'wnp'}
##'     calculates the within-node prevalence.
##' @param ... Additional arguments affecting the measure
##' @param i Indices specifying the nodes to include in the
##'     calculation of the prevalence. Default is \code{NULL}, which
##'     includes all nodes.
##' @param age For models with age categories, the age category to
##'     include in the calculation. Default is that all age categories
##'     are included.
##' @return Vector when type equals \code{'pop'} or \code{'bnp'} but
##'     matrix when type equals \code{'wnp'}.
##' @export
setGeneric("prevalence",
           function(model, type = c("pop", "bnp", "wnp"), i = NULL,
                    ...) standardGeneric("prevalence"))

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

##' Extract the events from a \code{SimInf_model} object
##'
##' @rdname events-methods
##' @param model The \code{model} to extract the events from.
##' @return \code{SimInf_events} object.
##' @export
##' @examples
##' ## Create an SIR model that includes scheduled events.
##' model <- SIR(u0     = u0_SIR(),
##'              tspan  = 1:(4 * 365),
##'              events = events_SIR(),
##'              beta   = 0.16,
##'              gamma  = 0.077)
##'
##' ## Extract the scheduled events from the model and
##' ## display summary
##' summary(events(model))
##'
##' ## Extract the scheduled events from the model and
##' ## plot summary
##' plot(events(model))
setGeneric("events", function(model) standardGeneric("events"))
