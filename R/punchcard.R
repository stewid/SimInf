## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015 - 2019  Stefan Engblom
## Copyright (C) 2015 - 2019  Stefan Widgren
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
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

##' Set a template for where to record result during a simualtion
##'
##' Using a sparse result matrix can save a lot of memory if the model
##' contains many nodes and time-points, but where only a few of the
##' data points are of interest for post-processing.
##'
##' Using a sparse result matrix can save a lot of memory if the model
##' contains many nodes and time-points, but where only a few of the
##' data points are of interest for post-processing. To use this
##' feature, a template has to be defined for which data points to
##' record. This is done using a \code{data.frame} that specifies the
##' time-points (column \sQuote{time}) and nodes (column
##' \sQuote{node}) to record the state of the compartments, see
##' \sQuote{Examples}. The specified time-points, nodes and
##' compartments must exist in the model, or an error is raised. Note
##' that specifying a template only affects which data-points are
##' recorded for post-processing, it does not affect how the solver
##' simulates the trajectory.
##' @param model The \code{model} to set a template for where to
##'     record result.
##' @param value A \code{data.frame} that specify the nodes,
##'     time-points and compartments to record the number of
##'     individuals at \code{tspan}. Use \code{NULL} to reset the
##'     model to record the number of inidividuals in each compartment
##'     in every node at each time-point in tspan.
##' @include check_arguments.R
##' @include SimInf_model.R
##' @export
##' @importFrom methods as
##' @importFrom Matrix sparseMatrix
##' @examples
##' ## Create an 'SIR' model with 6 nodes and initialize it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Run the model.
##' result <- run(model, threads = 1)
##'
##' ## Display the trajectory with data for every node at each
##' ## time-point in tspan.
##' trajectory(result)
##'
##' ## Assume we are only interested in nodes '2' and '4' at the
##' ## time-points '3' and '5'
##' df <- data.frame(time = c(3, 5, 3, 5),
##'                  node = c(2, 2, 4, 4),
##'                  S = c(TRUE, TRUE, TRUE, TRUE),
##'                  I = c(TRUE, TRUE, TRUE, TRUE),
##'                  R = c(TRUE, TRUE, TRUE, TRUE))
##' punchcard(model) <- df
##' result <- run(model, threads = 1)
##' trajectory(result)
##'
##' ## We can also specify to record only some of the compartments in
##' ## each time-step.
##' df <- data.frame(time = c(3, 5, 3, 5),
##'                  node = c(2, 2, 4, 4),
##'                  S = c(FALSE, TRUE, TRUE, TRUE),
##'                  I = c(TRUE, FALSE, TRUE, FALSE),
##'                  R = c(TRUE, FALSE, TRUE, TRUE))
##' punchcard(model) <- df
##' result <- run(model, threads = 1)
##' trajectory(result)
##'
##' ## A shortcut to specify to record all of the compartments in
##' ## each time-step is to only inlude node and time.
##' df <- data.frame(time = c(3, 5, 3, 5),
##'                  node = c(2, 2, 4, 4))
##' punchcard(model) <- df
##' result <- run(model, threads = 1)
##' trajectory(result)
##'
##' ## It is possible to use an empty 'data.frame' to specify
##' ## that no data-points should be recorded for the trajectory.
##' punchcard(model) <- data.frame()
##' result <- run(model, threads = 1)
##' trajectory(result)
##'
##' ## Use 'NULL' to reset the model to record data for every node at
##' ## each time-point in tspan.
##' punchcard(model) <- NULL
##' result <- run(model, threads = 1)
##' trajectory(result)
"punchcard<-" <- function(model, value)
{
    check_model_argument(model)

    dimU <- c(Nn(model) * Nc(model), length(model@tspan))
    dimV <- c(Nn(model) * Nd(model), length(model@tspan))
    iU <- numeric(0)
    jU <- numeric(0)
    iV <- numeric(0)
    jV <- numeric(0)
    x <- NA_integer_

    if (is.null(value)) {
        ## Create sparse templates.
        model@U_sparse <- sparseMatrix(i = iU, j = jU, x = x, dims = c(0, 0))
        model@V_sparse <- sparseMatrix(i = iV, j = jV, x = x, dims = c(0, 0))
        validObject(model)
        return(model)
    }

    if (!is.data.frame(value))
        stop("'value' argument is not a 'data.frame'")

    ## Clear dense result matrices.
    model@U <- matrix(data = integer(0), nrow = 0, ncol = 0)
    model@V <- matrix(data = numeric(0), nrow = 0, ncol = 0)

    if (nrow(value) == 0) {
        ## Create sparse templates.
        model@U_sparse <- sparseMatrix(i = iU, j = jU, x = x, dims = dimU)
        model@V_sparse <- sparseMatrix(i = iV, j = jV, x = x, dims = dimV)
        validObject(model)
        return(model)
    }

    ## Check the content in 'value'
    if (!all(c("node", "time") %in% names(value)))
        stop("'value' must have the columns 'time' and 'node'.")
    if (!is.numeric(value$time))
        value$time <- as.character(value$time)
    if (is.character(value$time))
        value$time <- model@tspan[match(value$time, names(model@tspan))]

    ## Sort the data.frame by time and node.
    value <- value[order(value$time, value$node), ]

    ## Match the nodes and time-points with the model.
    i <- match(value$node, seq_len(Nn(model)))
    if (any(is.na(i)))
        stop("Unable to match all nodes")
    j <- match(value$time, model@tspan)
    if (any(is.na(j)))
        stop("Unable to match all time-points to tspan")

    compartments <- setdiff(colnames(value), c("time", "node"))
    if (length(compartments) == 0) {
        ## Only node and time specified, add all compartments and mark
        ## them as TRUE.
        compartments <- c(rownames(model@S), rownames(model@v0))
        value[, compartments] <- TRUE
    }

    ## Coerce the compartments part of the data.frame to a logical
    ## vector that match the rows of compartments in the U matrix.
    if (any(compartments %in% rownames(model@S))) {
        valueU <- value[, c("time", "node", rownames(model@S))]
        valueU <- as.logical(t(as.matrix(valueU[, -(1:2)])))
        valueU[is.na(valueU)] <- FALSE

        ## Create an index to all of its compartments in the U
        ## matrix. Keep only compartments and time-points that are
        ## marked with TRUE.
        iU <- rep((i - 1) * Nc(model), each = Nc(model)) + seq_len(Nc(model))
        jU <- rep(j, each = Nc(model))
        iU <- iU[valueU]
        jU <- jU[valueU]
    }

    ## Coerce the compartments part of the data.frame to a logical
    ## vector that match the rows of continuous state compartments in
    ## the V matrix.
    if (any(compartments %in% rownames(model@v0))) {
        valueV <- value[, c("time", "node", rownames(model@v0))]
        valueV <- as.logical(t(as.matrix(valueV[, -(1:2)])))
        valueV[is.na(valueV)] <- FALSE

        ## Create an index to all of its compartments in the V
        ## matrix. Keep only compartments and time-points that are
        ## marked with TRUE.
        iV <- rep((i - 1) * Nd(model), each = Nd(model)) + seq_len(Nd(model))
        jV <- rep(j, each = Nd(model))
        iV <- iV[valueV]
        jV <- jV[valueV]
    }

    ## Create sparse templates.
    model@U_sparse <- sparseMatrix(i = iU, j = jU, x = x, dims = dimU)
    model@V_sparse <- sparseMatrix(i = iV, j = jV, x = x, dims = dimV)
    validObject(model)
    model
}
