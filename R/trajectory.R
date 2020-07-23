## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2020 Stefan Widgren
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

##' Determine if the trajectory is empty.
##' @noRd
do_is_trajectory_empty <- function(model, slots) {
    ## First, check the dimensions of each slot to determine if the
    ## trajectory is empty.
    empty <- all(vapply(slots, function(name) {
        all(identical(dim(slot(model, name)), c(0L, 0L)),
            identical(dim(slot(model, paste0(name, "_sparse"))), c(0L, 0L)))
    }, logical(1)))

    if (!isTRUE(empty)) {
        ## Need to also check the sparse slot.
        empty <- any(vapply(slots, function(name) {
            any(is.na(slot(model, paste0(name, "_sparse"))))
        }, logical(1)))
    }

    empty
}

##' Determine if the trajectory is empty.
##' @noRd
setGeneric(
    "is_trajectory_empty",
    signature = "model",
    function(model)
        standardGeneric("is_trajectory_empty"))

##' @include SimInf_model.R
##' @noRd
setMethod(
    "is_trajectory_empty",
    signature(model = "SimInf_model"),
    function(model) {
        do_is_trajectory_empty(model, c("U", "V"))
    }
)

##' Extract data in the internal matrix format
##'
##' @param m simulated data to extract.
##' @param n number of available compartments in the simulated data.
##' @param selected_compartments indices to selected compartments to
##'     extract from the simulated data and include in the matrix.
##' @param i subset of nodes to extract data from. If NULL, all
##'     available nodes are included.
##' @noRd
trajectory_as_is <- function(m, n, selected_compartments, i) {
    if (is.null(i)) {
        if (length(selected_compartments) == n)
            return(m)
        i <- seq_len(nrow(m) %/% n)
    }

    ## Extract subset of data.
    selected_compartments <- sort(selected_compartments)
    i <- rep(selected_compartments, length(i)) +
        rep((i - 1) * n, each = length(selected_compartments))
    m[i, seq_len(ncol(m)), drop = FALSE]
}

trajectory_data <- function(model, name) {
    x <- slot(model, paste0(name, "_sparse"))
    if (!identical(dim(x), c(0L, 0L)))
        return(x)
    slot(model, name)
}

##' Extract data from a simulated trajectory
##'
##' Extract the number of individuals in each compartment in every
##' node after generating a single stochastic trajectory with
##' \code{\link{run}}.
setGeneric(
    "trajectory",
    signature = "model",
    function(model, ...)
        standardGeneric("trajectory"))

##' @rdname trajectory
##' @section Internal format of the discrete state variables:
##'     Description of the layout of the internal matrix (\code{U})
##'     that is returned if \code{as.is = TRUE}. \code{U[, j]}
##'     contains the number of individuals in each compartment at
##'     \code{tspan[j]}. \code{U[1:Nc, j]} contains the number of
##'     individuals in node 1 at \code{tspan[j]}. \code{U[(Nc + 1):(2
##'     * Nc), j]} contains the number of individuals in node 2 at
##'     \code{tspan[j]} etc, where \code{Nc} is the number of
##'     compartments in the model. The dimension of the matrix is
##'     \eqn{N_n N_c \times} \code{length(tspan)} where \eqn{N_n} is
##'     the number of nodes.
##' @section Internal format of the continuous state variables:
##'     Description of the layout of the matrix that is returned if
##'     \code{as.is = TRUE}. The result matrix for the real-valued
##'     continuous state. \code{V[, j]} contains the real-valued state
##'     of the system at \code{tspan[j]}. The dimension of the matrix
##'     is \eqn{N_n}\code{dim(ldata)[1]} \eqn{\times}
##'     \code{length(tspan)}.
##' @param model the \code{model} to extract the result from.
##' @param compartments specify the names of the compartments to
##'     extract data from. The compartments can be specified as a
##'     character vector e.g. \code{compartments = c('S', 'I', 'R')},
##'     or as a formula e.g. \code{compartments = ~S+I+R} (see
##'     \sQuote{Examples}). Default (\code{compartments=NULL}) is to
##'     extract the number of individuals in each compartment i.e. the
##'     data from all discrete state compartments in the model. In
##'     models that also have continuous state variables e.g. the
##'     \code{SISe} model, use \code{~.} instead of \code{NULL} to
##'     also include these.
##' @param node indices specifying the subset of nodes to include when
##'     extracting data. Default (\code{node = NULL}) is to extract data
##'     from all nodes.
##' @param as.is the default (\code{as.is = FALSE}) is to generate a
##'     \code{data.frame} with one row per node and time-step with the
##'     number of individuals in each compartment. Using \code{as.is =
##'     TRUE} returns the result as a matrix, which is the internal
##'     format (see \sQuote{Details}).
##' @param ... Additional arguments. Not used.
##' @return A \code{data.frame} if \code{as.is = FALSE}, else a
##'     matrix.
##' @include SimInf_model.R
##' @include check_arguments.R
##' @include match_compartments.R
##' @include prevalence.R
##' @export
##' @importFrom methods is
##' @examples
##' ## Create an 'SIR' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##'
##' ## Extract the number of individuals in each compartment at the
##' ## time-points in 'tspan'.
##' trajectory(result)
##'
##' ## Extract the number of recovered individuals in the first node
##' ## at the time-points in 'tspan'.
##' trajectory(result, compartments = "R", node = 1)
##'
##' ## Extract the number of recovered individuals in the first and
##' ## third node at the time-points in 'tspan'.
##' trajectory(result, compartments = "R", node = c(1, 3))
##'
##' ## Create an 'SISe' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6)
##' model <- SISe(u0 = u0, tspan = 1:10, phi = rep(0, 6),
##'     upsilon = 0.02, gamma = 0.1, alpha = 1, epsilon = 1.1e-5,
##'     beta_t1 = 0.15, beta_t2 = 0.15, beta_t3 = 0.15, beta_t4 = 0.15,
##'     end_t1 = 91, end_t2 = 182, end_t3 = 273, end_t4 = 365)
##'
##' ## Run the model
##' result <- run(model)
##'
##' ## Extract the continuous state variable 'phi' which represents
##' ## the environmental infectious pressure.
##' trajectory(result, "phi")
setMethod(
    "trajectory",
    signature(model = "SimInf_model"),
    function(model, compartments = NULL, node = NULL, as.is = FALSE, ...) {
        if (is_trajectory_empty(model)) {
            stop("Please run the model first, the trajectory is empty.",
                 call. = FALSE)
        }

        compartments <- match_compartments(compartments = compartments,
                                           ok_combine = !isTRUE(as.is),
                                           ok_lhs = FALSE,
                                           U = rownames(model@S),
                                           V = rownames(model@v0))

        node <- check_node_argument(model, node)

        if (isTRUE(as.is)) {
            ## Extract data in the internal matrix format.
            if (length(compartments$rhs$U)) {
                return(trajectory_as_is(trajectory_data(model, "U"), Nc(model),
                                        compartments$rhs$U, node))
            }

            return(trajectory_as_is(trajectory_data(model, "V"), Nd(model),
                                    compartments$rhs$V, node))
        }

        ## Coerce the dense/sparse 'U' and 'V' matrices to a
        ## data.frame with one row per node and time-point with data
        ## from the specified discrete and continuous states.
        .Call(SimInf_trajectory,
              trajectory_data(model, "U"), compartments$rhs$U,
              attr(compartments$rhs$U, "available_compartments"),
              trajectory_data(model, "V"), compartments$rhs$V,
              attr(compartments$rhs$V, "available_compartments"),
              model@tspan, n_nodes(model), node, "node")
    }
)
