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

        ## Cannot combine data from U and V when as.is = TRUE or when
        ## both U and V are sparse.
        if (all(!is.null(compartments_U), !is.null(compartments_V))) {
            if (isTRUE(as.is))
                stop("Select either continuous or discrete compartments")
            if (all(!identical(dim(model@U_sparse), c(0L, 0L)),
                    !identical(dim(model@V_sparse), c(0L, 0L))))
                stop("Select either continuous or discrete compartments")
        }
    }

    ## Check the 'node' argument
    node <- check_node_argument(model, node)

    ## The arguments seem ok...go on and extract the trajectory

    ## Check to extract sparse data from V
    if (!identical(dim(model@V_sparse), c(0L, 0L))) {
        if (!is.null(compartments_V)) {
            if (isTRUE(as.is))
                return(model@V_sparse)

            ## Coerce the sparse 'V_sparse' matrix to a data.frame with
            ## one row per node and time-point with the values of the
            ## continuous state variables.
            return(sparse2df(model@V_sparse, Nd(model), model@tspan,
                             rownames(model@v0), NA_real_))
        }
    }

    ## Check to extract sparse data from U
    if (!identical(dim(model@U_sparse), c(0L, 0L))) {
        if (isTRUE(as.is))
            return(model@U_sparse)

        ## Coerce the sparse 'U_sparse' matrix to a data.frame with
        ## one row per node and time-point with the number of
        ## individuals in each compartment.
        return(sparse2df(model@U_sparse, Nc(model),
                         model@tspan, rownames(model@S)))
    }

    ## Check to extract data in internal matrix format
    if (isTRUE(as.is)) {
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
            j <- rep(compartments_V, length(node))
            j <- j + rep((node - 1) * Nd(model), each = length(compartments_V))
            return(model@V[j, , drop = FALSE])
        }

        ## Extract subset of data from U
        compartments_U <- sort(match(compartments_U, rownames(model@S)))
        j <- rep(compartments_U, length(node))
        j <- j + rep((node - 1) * Nc(model), each = length(compartments_U))
        return(model@U[j, , drop = FALSE])
    }

    ## Coerce the dense 'U' and 'V' matrices to a data.frame with one
    ## row per node and time-point with data from the specified
    ## discrete and continuous states.
    mU <- NULL
    mV <- NULL

    ## Handle first cases where all data in U and/or V are extracted,
    ## where 'all' indicates that all compartments in U or V are
    ## specified.
    ##
    ## compartments_U compartments_V output
    ##     NULL           NULL         U
    ##     NULL            all         V
    ##      all           NULL         U
    ##      all            all        U+V
    if (is.null(node)) {
        if (is.null(compartments_U)) {
            if (is.null(compartments_V)) {
                mU <- matrix(as.integer(model@U), ncol = Nc(model), byrow = TRUE)
            } else if (identical(length(compartments_V), Nd(model))) {
                mV <- matrix(as.numeric(model@V), ncol = Nd(model), byrow = TRUE)
            }
        } else if (identical(length(compartments_U), Nc(model))) {
            if (is.null(compartments_V)) {
                mU <- matrix(as.integer(model@U), ncol = Nc(model), byrow = TRUE)
            } else if (identical(length(compartments_V), Nd(model))) {
                mU <- matrix(as.integer(model@U), ncol = Nc(model), byrow = TRUE)
                mV <- matrix(as.numeric(model@V), ncol = Nd(model), byrow = TRUE)
            }
        }

        if (!is.null(mU))
            colnames(mU) <- rownames(model@S)
        if (!is.null(mV))
            colnames(mV) <- rownames(model@v0)
    }

    ## Handle cases where a subset of data in U and/or V are
    ## extracted.
    if (all(is.null(mU), is.null(mV))) {
        if (is.null(node))
            node <- seq_len(Nn(model))

        if (all(is.null(compartments_U), is.null(compartments_V)))
            compartments_U <- rownames(model@S)

        if (!is.null(compartments_U)) {
            ## Extract a subset of data from U
            compartments_U <- sort(match(compartments_U, rownames(model@S)))
            j <- rep(compartments_U, length(node))
            j <- j + rep((node - 1) * Nc(model), each = length(compartments_U))
            k <- (seq_len(length(model@tspan)) - 1) * Nc(model) * Nn(model)
            k <- rep(k, each = length(j))
            j <- rep(j, length(model@tspan))
            j <- j + k
            mU <- matrix(as.integer(model@U[j]),
                         ncol = length(compartments_U),
                         byrow = TRUE)
            colnames(mU) <- rownames(model@S)[compartments_U]
        }

        if (!is.null(compartments_V)) {
            ## Extract a subset of data from V
            compartments_V <- sort(match(compartments_V, rownames(model@v0)))
            j <- rep(compartments_V, length(node))
            j <- j + rep((node - 1) * Nd(model), each = length(compartments_V))
            k <- (seq_len(length(model@tspan)) - 1) * Nd(model) * Nn(model)
            k <- rep(k, each = length(j))
            j <- rep(j, length(model@tspan))
            j <- j + k
            mV <- matrix(as.numeric(model@V[j]),
                         ncol = length(compartments_V),
                         byrow = TRUE)
            colnames(mV) <- rownames(model@v0)[compartments_V]
        }
    }

    if (is.null(node))
        node = seq_len(Nn(model))

    time <- names(model@tspan)
    if (is.null(time))
        time <- as.integer(model@tspan)
    time <- rep(time, each = length(node))

    result <- data.frame(node = node, time = time, stringsAsFactors = FALSE)
    if (!is.null(mU))
        result <- cbind(result, as.data.frame(mU))

    if (!is.null(mV))
        result <- cbind(result, as.data.frame(mV))

    result
}
