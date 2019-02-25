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

##' Coerce a sparse matrix to a data.frame
##'
##' Utility function to coerce a sparse matrix (U_sparse or V_sparse)
##' to a data.frame.
##' @param m sparse matrix to coerce.
##' @param n number of rows per node.
##' @param tspan time points in trajectory.
##' @param lbl labels for data e.g. compartments.
##' @param value default value.
##' @return \code{data.frame}
##' @noRd
sparse2df <- function(m, n, tspan, lbl, value = NA_integer_) {
    ## Determine nodes and time-points with output.
    node <- as.integer(ceiling((m@i + 1) / n))
    time <- names(tspan)
    if (is.null(time))
        time <- as.integer(tspan)
    time <- cbind(time, diff(m@p))
    time <- unlist(apply(time, 1, function(x) rep(x[1], x[2])))

    ## Determine unique combinations of node and time
    i <- !duplicated(cbind(node, time))
    node <- node[i]
    time <- time[i]

    ## Use node and time to determine the required size
    ## of a matrix to hold all output data and fill it
    ## with NA values.
    values <- matrix(value, nrow = sum(i), ncol = n)
    colnames(values) <- lbl

    ## And then update non-NA items with values from m.
    i <- cumsum(i)
    j <- m@i %% n + 1
    if (is.integer(value)) {
        values[matrix(c(i, j), ncol = 2)] <- as.integer(m@x)
    } else {
        values[matrix(c(i, j), ncol = 2)] <- m@x
    }

    cbind(node = node,
          time = time,
          as.data.frame(values),
          stringsAsFactors = FALSE)
}

denseU2df <- function(model, compartments, node)
{
    mU <- NULL
    if (is.null(node)) {
        if (length(compartments) == Nc(model)) {
            mU <- matrix(as.integer(model@U), ncol = Nc(model), byrow = TRUE)
            colnames(mU) <- rownames(model@S)
        }

        node <- seq_len(Nn(model))
    }

    if (is.null(mU)) {
        ## Extract a subset of data from U
        compartments <- sort(match(compartments, rownames(model@S)))
        i <- rep(compartments, length(node))
        i <- i + rep((node - 1) * Nc(model), each = length(compartments))
        j <- (seq_len(length(model@tspan)) - 1) * Nc(model) * Nn(model)
        j <- rep(j, each = length(i))
        i <- rep(i, length(model@tspan))
        i <- i + j
        mU <- matrix(as.integer(model@U[i]),
                     ncol = length(compartments),
                     byrow = TRUE)
        colnames(mU) <- rownames(model@S)[compartments]
    }

    time <- names(model@tspan)
    if (is.null(time))
        time <- as.integer(model@tspan)
    time <- rep(time, each = length(node))

    cbind(data.frame(node = node, time = time, stringsAsFactors = FALSE),
          as.data.frame(mU))
}

denseV2df <- function(model, compartments, node)
{
    mV <- NULL
    if (is.null(node)) {
        if (length(compartments) == Nd(model)) {
            mV <- matrix(as.numeric(model@V), ncol = Nd(model), byrow = TRUE)
            colnames(mV) <- rownames(model@v0)
        }

        node <- seq_len(Nn(model))
    }

    if (is.null(mV)) {
        ## Extract a subset of data from V
        compartments <- sort(match(compartments, rownames(model@v0)))
        i <- rep(compartments, length(node))
        i <- i + rep((node - 1) * Nd(model), each = length(compartments))
        j <- (seq_len(length(model@tspan)) - 1) * Nd(model) * Nn(model)
        j <- rep(j, each = length(i))
        i <- rep(i, length(model@tspan))
        i <- i + j
        mV <- matrix(as.numeric(model@V[i]),
                     ncol = length(compartments),
                     byrow = TRUE)
        colnames(mV) <- rownames(model@v0)[compartments]
    }

    time <- names(model@tspan)
    if (is.null(time))
        time <- as.integer(model@tspan)
    time <- rep(time, each = length(node))

    cbind(data.frame(node = node, time = time, stringsAsFactors = FALSE),
          as.data.frame(mV))
}

##' Determine if the trajectory is empty.
##' @noRd
is_trajectory_empty <- function(model)
{
    if (all(identical(dim(model@U), c(0L, 0L)),
            identical(dim(model@U_sparse), c(0L, 0L)),
            identical(dim(model@V), c(0L, 0L)),
            identical(dim(model@V_sparse), c(0L, 0L)))) {
        return(TRUE)
    }

    if (any(is.na(model@U_sparse@x)) || any(is.na(model@V_sparse@x)))
        return(TRUE)

    FALSE
}

##' Extract data in internal matrix format.
##' @noRd
trajectory_as_is <- function(model, compartments_U, compartments_V, node)
{
    if (!identical(dim(model@V_sparse), c(0L, 0L))) {
        if (!is.null(compartments_V))
            return(model@V_sparse)
    }

    if (!identical(dim(model@U_sparse), c(0L, 0L)))
        return(model@U_sparse)

    if (is.null(node)) {
        if (is.null(compartments_U)) {
            if (is.null(compartments_V))
                return(model@U)
            if (identical(length(compartments_V), Nd(model)))
                return(model@V)
        } else if (identical(length(compartments_U), Nc(model))) {
            return(model@U)
        }
    }

    if (is.null(node))
        node <- seq_len(Nn(model))

    if (all(is.null(compartments_U), is.null(compartments_V)))
        compartments_U <- rownames(model@S)

    if (is.null(compartments_U)) {
        ## Extract subset of data from V
        compartments_V <- sort(match(compartments_V, rownames(model@v0)))
        i <- rep(compartments_V, length(node))
        i <- i + rep((node - 1) * Nd(model), each = length(compartments_V))
        return(model@V[i, , drop = FALSE])
    }

    ## Extract subset of data from U
    compartments_U <- sort(match(compartments_U, rownames(model@S)))
    i <- rep(compartments_U, length(node))
    i <- i + rep((node - 1) * Nc(model), each = length(compartments_U))
    model@U[i, , drop = FALSE]
}

##' Extract data from a simulated trajectory
##'
##' Extract the number of individuals in each compartment in every
##' node after generating a single stochastic trajectory with
##' \code{\link{run}}.
##'
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
##' @return A \code{data.frame} if \code{as.is = FALSE}, else a
##'     matrix.
##' @include SimInf_model.R
##' @include check_arguments.R
##' @include prevalence.R
##' @export
##' @importFrom methods is
##' @importFrom stats terms
##' @examples
##' ## Create an 'SIR' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model, threads = 1)
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
##' result <- run(model, threads = 1)
##'
##' ## Extract the continuous state variable 'phi' which represents
##' ## the environmental infectious pressure.
##' trajectory(result, "phi")
trajectory <- function(model, compartments = NULL, node = NULL, as.is = FALSE)
{
    ## Check that the arguments are ok...

    check_model_argument(model)

    if (is_trajectory_empty(model))
        stop("Please run the model first, the trajectory is empty")

    ## Split the 'compartments' argument to match the compartments in
    ## U and V.
    if (is(compartments, "formula")) {
        compartments <- as.character(compartments)
        if (!identical(length(compartments), 2L))
            stop("Invalid formula specification of 'compartments'")
        compartments <- parse_formula_item(
            compartments[2], c(rownames(model@S), rownames(model@v0)))
    }

    compartments_U <- NULL
    compartments_V <- NULL
    if (!is.null(compartments)) {
        compartments <- unique(as.character(compartments))

        ## Match compartments in U
        lbl <- rownames(model@S)
        compartments_U <- compartments[compartments %in% lbl]
        if (length(compartments_U) > 0) {
            compartments <- setdiff(compartments, compartments_U)
        } else {
            compartments_U <- NULL
        }

        ## Match compartments in V
        if (length(compartments) > 0) {
            if (Nd(model) > 0) {
                lbl <- rownames(model@v0)
                compartments_V <- compartments[compartments %in% lbl]
                if (length(compartments_V) > 0) {
                    compartments <- setdiff(compartments, compartments_V)
                } else {
                    compartments_V <- NULL
                }
            }
        }

        if (length(compartments) > 0) {
            stop("Non-existing compartment(s) in model: ",
                 paste0("'", compartments, "'", collapse = ", "))
        }

        ## Cannot combine data from U and V when as.is = TRUE.
        if (!is.null(compartments_U) && !is.null(compartments_V) && isTRUE(as.is))
            stop("Select either continuous or discrete compartments")
    }

    if (is.null(compartments_U) && is.null(compartments_V))
        compartments_U <- rownames(model@S)

    ## Check the 'node' argument
    node <- check_node_argument(model, node)

    ## The arguments seem ok...go on and extract the trajectory

    ## Check to extract data in internal matrix format
    if (isTRUE(as.is))
        return(trajectory_as_is(model, compartments_U, compartments_V, node))

    ## Coerce the dense/sparse 'U' and 'V' matrices to a data.frame
    ## with one row per node and time-point with data from the
    ## specified discrete and continuous states.

    dfV <- NULL
    if (!is.null(compartments_V)) {
        if (identical(dim(model@V_sparse), c(0L, 0L))) {
            dfV <- denseV2df(model, compartments_V, node)
        } else {
            dfV <- sparse2df(model@V_sparse, Nd(model), model@tspan,
                             rownames(model@v0), NA_real_)
        }
    }

    dfU <- NULL
    if (!is.null(compartments_U)) {
        if (identical(dim(model@U_sparse), c(0L, 0L))) {
            dfU <- denseU2df(model, compartments_U, node)
        } else {
            dfU <- sparse2df(model@U_sparse, Nc(model),
                             model@tspan, rownames(model@S))
        }
    }

    if (is.null(dfV))
        return(dfU)
    if (is.null(dfU))
        return(dfV)
    merge(dfU, dfV, all = TRUE, sort = FALSE)
}
