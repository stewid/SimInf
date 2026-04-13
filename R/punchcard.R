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

##' Set a sparse recording template for simulation results
##'
##' Define a template (or "punchcard") specifying which data points
##' (nodes, time-points, and compartments) should be recorded during a
##' simulation. This feature is useful for saving memory when the
##' model has many nodes or time-points, but only a subset of the data
##' is needed for post-processing.
##'
##' The template is specified as a \code{data.frame} with columns:
##' \itemize{
##'   \item \code{time}: The time-points to record.
##'   \item \code{node}: The node indices to record.
##'   \item \code{CompartmentName}: Logical columns (e.g., \code{S},
##'     \code{I}, \code{R}) indicating whether to record that specific
##'     compartment for the corresponding row. If a compartment column
##'     is omitted, it defaults to \code{FALSE} (not recorded). If
##'     only \code{time} and \code{node} are provided, all
##'     compartments are recorded for those points.
##' }
##'
##' \strong{Important:} The specified \code{time} values, \code{node}
##' indices, and compartment names must exist in the model. If any
##' value in the template does not match the model's configuration
##' (e.g., a time point not in \code{tspan} or a compartment not
##' defined in the model), an \strong{error will be raised} when the
##' template is set.
##'
##' Note that the template only affects which data is stored in the
##' result object; it does not affect how the solver simulates the
##' trajectory.
##'
##' @param model A \code{SimInf_model} object.
##' @param value A \code{data.frame} defining the recording template,
##'     or \code{NULL} to reset the model to record all compartments
##'     for all nodes at all time-points in \code{tspan}.
##' @include check_arguments.R
##' @include SimInf_model.R
##' @export
##' @examples
##' ## For reproducibility, set the seed and number of threads.
##' set.seed(123)
##' set_num_threads(1)
##'
##' ## Create an 'SIR' model with 6 nodes
##' u0 <- data.frame(
##'   S = 100:105,
##'   I = 1:6,
##'   R = rep(0, 6)
##' )
##' model <- SIR(
##'   u0 = u0,
##'   tspan = 1:10,
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Run the model with default recording (all data)
##' result_full <- run(model)
##' head(trajectory(result_full))
##'
##' ## Define a template to record only nodes 2 and 4 at times 3 and 5
##' ## for all compartments.
##' df <- data.frame(
##'   time = c(3, 5, 3, 5),
##'   node = c(2, 2, 4, 4),
##'   S = TRUE, I = TRUE, R = TRUE
##' )
##' punchcard(model) <- df
##'
##' result_sparse <- run(model)
##' trajectory(result_sparse)
##'
##' ## Record only specific compartments (e.g., S and R, but not I)
##' df <- data.frame(
##'   time = c(3, 5, 3, 5),
##'   node = c(2, 2, 4, 4),
##'   S = TRUE, I = FALSE, R = TRUE
##' )
##' punchcard(model) <- df
##' result_partial <- run(model)
##' trajectory(result_partial)
##'
##' ## Shortcut: If only 'time' and 'node' are specified, all
##' ## compartments are recorded.
##' df <- data.frame(
##'   time = c(3, 5, 3, 5),
##'   node = c(2, 2, 4, 4)
##' )
##' punchcard(model) <- df
##'
##' ## Reset to record all data (equivalent to no template)
##' punchcard(model) <- NULL
##' result_reset <- run(model)
##' head(trajectory(result_reset))
## nolint start: brace_linter
setGeneric(
    "punchcard<-",
    signature = "model",
    function(model,
             value)
        standardGeneric("punchcard<-")
)
## nolint end

##' @rdname punchcard-set
##' @export
setMethod(
    "punchcard<-",
    signature(model = "SimInf_model"),
    function(model, value) {
        template <- create_template(value, model@tspan, seq_len(n_nodes(model)),
                                    rownames(model@S), integer(0))
        model@U <- template$dense
        model@U_sparse <- template$sparse

        template <- create_template(value, model@tspan, seq_len(n_nodes(model)),
                                    rownames(model@v0), numeric(0))
        model@V <- template$dense
        model@V_sparse <- template$sparse

        methods::validObject(model)
        model
    }
)

##' Create  template for where to record result during a simualtion
##'
##' @param value A \code{data.frame} that specify the nodes,
##'     time-points and compartments to record the number of
##'     individuals at \code{tspan}. Use \code{NULL} to reset the
##'     model to record the number of inidividuals in each compartment
##'     in every node at each time-point in tspan.
##' @param tspan time points in trajectory.
##' @param compartments available compartments in the simulated data.
##' @param data default data in dense matrix.
##' @noRd
create_template <- function(value, tspan, nodes, compartments, data) {
    if (is.null(value)) {
        dense <- matrix(data = data, nrow = 0, ncol = 0)
        sparse <- methods::new("dgCMatrix")
        return(list(dense = dense, sparse = sparse))
    }

    if (!is.data.frame(value))
        stop("'value' argument is not a 'data.frame'.", call. = FALSE)

    if (nrow(value) == 0) {
        dense <- matrix(data = data, nrow = 0, ncol = 0)
        dims <- c(length(nodes) * length(compartments), length(tspan))
        sparse <- Matrix::sparseMatrix(i = numeric(0), j = numeric(0),
                                       x = NA_real_, dims = dims)
        return(list(dense = dense, sparse = sparse))
    }

    ## Check the content in 'value'
    if (!all(c("node", "time") %in% names(value)))
        stop("'value' must have the columns 'time' and 'node'.", call. = FALSE)
    if (!is.numeric(value$time))
        value$time <- as.character(value$time)
    if (is.character(value$time))
        value$time <- tspan[match(value$time, names(tspan))]

    ## Sort the data.frame by time and node.
    value <- value[order(value$time, value$node), ]

    ## Match the nodes and time-points with the model.
    i <- match(value$node, nodes)
    if (anyNA(i))
        stop("Unable to match all nodes.", call. = FALSE)
    j <- match(value$time, tspan)
    if (anyNA(j))
        stop("Unable to match all time-points to tspan.", call. = FALSE)

    selected_compartments <- setdiff(colnames(value), c("time", "node"))
    if (length(selected_compartments) == 0) {
        ## Only node and time specified, select all compartments and
        ## mark them as TRUE.
        selected_compartments <- compartments
        value[, selected_compartments] <- TRUE
    }

    ## Coerce the compartments part of the data.frame to a logical
    ## vector that match the rows of compartments in the matrix.
    if (any(selected_compartments %in% compartments)) {
        value <- value[, c("time", "node", compartments)]
        value <- as.logical(t(as.matrix(value[, c(-1, -2)])))
        value[is.na(value)] <- FALSE

        ## Create an index to all of its compartments in the
        ## matrix. Keep only compartments and time-points that are
        ## marked with TRUE.
        i <- rep((i - 1) * length(compartments),
                 each = length(compartments)) + seq_along(compartments)
        i <- i[value]
        j <- rep(j, each = length(compartments))
        j <- j[value]

        dims <- c(length(nodes) * length(compartments), length(tspan))
        d1_times_d2 <- as.numeric(dims[1]) * as.numeric(dims[2])
        if (sum(value, na.rm = TRUE) == d1_times_d2) {
            dense <- matrix(data = data, nrow = 0, ncol = 0)
            sparse <- methods::new("dgCMatrix")
        } else {
            dense <- matrix(data = data, nrow = 0, ncol = 0)
            sparse <- Matrix::sparseMatrix(i = i, j = j, x = NA_real_,
                                           dims = dims)
        }

        return(list(dense = dense, sparse = sparse))
    }

    dense <- matrix(data = data, nrow = 0, ncol = 0)
    dims <- c(length(nodes) * length(compartments), length(tspan))
    sparse <- Matrix::sparseMatrix(i = numeric(0), j = numeric(0),
                                   x = NA_real_, dims = dims)
    list(dense = dense, sparse = sparse)
}
