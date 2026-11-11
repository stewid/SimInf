## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2026 Stefan Widgren
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
            anyNA(slot(model, paste0(name, "_sparse")))
        }, logical(1)))
    }

    empty
}

##' Determine if the trajectory is empty.
##' @noRd
## nolint start: brace_linter
setGeneric(
    "is_trajectory_empty",
    signature = "model",
    function(model)
        standardGeneric("is_trajectory_empty")
)
## nolint end

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
##' @param index indices specifying the subset of nodes to include
##'     when extracting data. Default (\code{index = NULL}) is to
##'     extract data from all nodes.
##' @noRd
trajectory_as_is <- function(m, n, selected_compartments, index) {
    if (is.null(index)) {
        if (length(selected_compartments) == n)
            return(m)
        index <- seq_len(nrow(m) %/% n)
    }

    ## Extract subset of data.
    selected_compartments <- sort(selected_compartments)
    index <- rep(selected_compartments, length(index)) +
        rep((index - 1) * n, each = length(selected_compartments))
    m[index, seq_len(ncol(m)), drop = FALSE]
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
##' \code{\link{run}}. For a \code{\linkS4class{SimInf_pfilter}}
##' object, extract the filtered trajectory from running a particle
##' filter.
##' @section Internal format of the discrete state variables: When
##'     \code{format = "matrix"}, the result is the internal matrix
##'     \code{U}. The column \code{U[, j]} contains the number of
##'     individuals in each compartment at \code{tspan[j]}, where the
##'     first \eqn{N_c} rows contain node 1, the next \eqn{N_c} rows
##'     contain node 2, and so on, where \eqn{N_c} is the number of
##'     compartments in the model. Consequently, the dimension of the
##'     matrix is \eqn{N_n N_c \times} \code{length(tspan)}, where
##'     \eqn{N_n} is the number of nodes. The internal format also
##'     allows replicates of each node. Each replicate adds new
##'     columns to \code{U}, so that the data for each replicate is in
##'     blocks of \code{length(tspan)} columns.
##' @section Internal format of the continuous state variables: When
##'     \code{format = "matrix"}, the result is the internal matrix
##'     \code{V}. The column \code{V[, j]} contains the real-valued
##'     continuous state of the system at \code{tspan[j]}, where the
##'     rows are grouped per node and the number of rows per node is
##'     determined by the number of continuous state variables in the
##'     model. Consequently, the dimension of the matrix is \eqn{N_n
##'     N_d \times} \code{length(tspan)}, where \eqn{N_n} is the
##'     number of nodes and \eqn{N_d} is the number of continuous
##'     state variables. The internal format also allows replicates of
##'     each node. Each replicate adds new columns to \code{V}, so
##'     that the data for each replicate is in blocks of
##'     \code{length(tspan)} columns.
##' @param model the \code{SimInf_model} object to extract the result
##'     from.
##' @template compartments-param
##' @template index-param
##' @param format the default (\code{format = "data.frame"}) is to
##'     generate a \code{data.frame} with one row per node and
##'     time-step with the number of individuals in each
##'     compartment. When the model contains multiple replicates of
##'     each node, the \code{data.frame} also contains one column
##'     \code{replicate}. Using \code{format = "matrix"} returns the
##'     result as a matrix, which is the internal format (see
##'     \sQuote{Details}).
##' @return A \code{data.frame} if \code{format = "data.frame"}, else
##'     a matrix.
##' @seealso \code{\link[=run,SimInf_model-method]{run}} for
##'     generating a trajectory, and
##'     \code{\link[=prevalence,SimInf_model-method]{prevalence}} for
##'     calculating prevalence from trajectory data.  The vignette
##'     "Post-process data in a trajectory" for a tutorial on
##'     post-processing, see
##'     \code{vignette("Post-process data in a trajectory", package =
##'     "SimInf")}.
##' @examples
##' ## For reproducibility, specify the number of threads.
##' set_num_threads(1)
##'
##' ## Create an 'SIR' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' model <- SIR(
##'   u0 = data.frame(
##'     S = 100:105,
##'     I = 1:6,
##'     R = rep(0, 6)),
##'   tspan = 1:10,
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Run the model with a fixed seed for reproducibility.
##' result <- run(model, seed = 123)
##'
##' ## Extract the number of individuals in each compartment at the
##' ## time-points in 'tspan'.
##' trajectory(result)
##'
##' ## Extract the number of recovered individuals in the first node
##' ## at the time-points in 'tspan'.
##' trajectory(result, compartments = "R", index = 1)
##'
##' ## Extract the number of recovered individuals in the first and
##' ## third node at the time-points in 'tspan'.
##' trajectory(result, compartments = "R", index = c(1, 3))
##'
##' ## Create an 'SISe' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' model <- SISe(
##'   u0 = data.frame(
##'     S = 100:105,
##'     I = 1:6),
##'   tspan = 1:10,
##'   phi = rep(0, 6),
##'   upsilon = 0.02,
##'   gamma = 0.1,
##'   alpha = 1,
##'   epsilon = 1.1e-5,
##'   beta_t1 = 0.15,
##'   beta_t2 = 0.15,
##'   beta_t3 = 0.15,
##'   beta_t4 = 0.15,
##'   end_t1 = 91,
##'   end_t2 = 182,
##'   end_t3 = 273,
##'   end_t4 = 365
##' )
##'
##' ## Run the model with a fixed seed for reproducibility.
##' result <- run(model, seed = 123)
##'
##' ## Extract the continuous state variable 'phi' which represents
##' ## the environmental infectious pressure.
##' trajectory(result, compartments = "phi")
## nolint start: brace_linter
setGeneric(
    "trajectory",
    signature = "model",
    function(model,
             compartments = NULL,
             index = NULL,
             format = c("data.frame", "matrix"))
        standardGeneric("trajectory")
)
## nolint end

##' @rdname trajectory
##' @include SimInf_model.R
##' @include check_arguments.R
##' @include match_compartments.R
##' @export
setMethod(
    "trajectory",
    signature(model = "SimInf_model"),
    function(model, compartments, index, format) {
        if (is_trajectory_empty(model)) {
            stop("Please run the model first, the trajectory is empty.",
                 call. = FALSE)
        }

        format <- match.arg(format)
        compartments <- match_compartments(compartments = compartments,
                                           ok_combine =
                                               identical(format, "data.frame"),
                                           ok_lhs = FALSE,
                                           U = rownames(model@S),
                                           V = rownames(model@v0))

        index <- check_node_index_argument(model, index)

        if (identical(format, "matrix")) {
            ## Extract data in the internal matrix format.
            if (length(compartments$rhs$U)) {
                return(trajectory_as_is(trajectory_data(model, "U"), Nc(model),
                                        compartments$rhs$U, index))
            }

            return(trajectory_as_is(trajectory_data(model, "V"), Nd(model),
                                    compartments$rhs$V, index))
        }

        ## Coerce the dense/sparse 'U' and 'V' matrices to a
        ## data.frame with one row per node and time-point with data
        ## from the specified discrete and continuous states.
        .Call(SimInf_trajectory,
              trajectory_data(model, "U"), compartments$rhs$U,
              attr(compartments$rhs$U, "available_compartments"),
              trajectory_data(model, "V"), compartments$rhs$V,
              attr(compartments$rhs$V, "available_compartments"),
              model@tspan, n_nodes(model), index, "node",
              n_replicates(model))
    }
)

##' @rdname trajectory
##' @include pfilter.R
##' @export
setMethod(
    "trajectory",
    signature(model = "SimInf_pfilter"),
    function(model, compartments, index, format) {
        trajectory(model@model, compartments, index, format)
    }
)
