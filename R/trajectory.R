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

## Split the 'compartments' argument to match the compartments in U
## and V.
match_compartments <- function(model, compartments, as.is)
{
    U <- NULL
    V <- NULL

    if (!is.null(compartments)) {
        compartments <- unique(as.character(compartments))

        ## Match compartments in U
        i <- compartments %in% rownames(model@S)
        if (any(i))
            U <- compartments[i]

        ## Match compartments in V
        i <- compartments %in% rownames(model@v0)
        if (any(i))
            V <- compartments[i]

        compartments <- setdiff(compartments, c(U, V))
        if (length(compartments) > 0) {
            stop("Non-existing compartment(s) in model: ",
                 paste0("'", compartments, "'", collapse = ", "),
                 ".", call. = FALSE)
        }

        ## Cannot combine data from U and V when as.is = TRUE.
        if (!is.null(U) && !is.null(V) && isTRUE(as.is)) {
            stop("Select either continuous or discrete compartments.",
                 call. = FALSE)
        }
    }

    if (is.null(U) && is.null(V))
        U <- rownames(model@S)

    list(U = U, V = V)
}

parse_formula <- function(model, compartments)
{
    compartments <- as.character(compartments)
    if (!identical(length(compartments), 2L))
        stop("Invalid formula specification of 'compartments'.", call. = FALSE)

    parse_formula_item(
        compartments[2],
        c(rownames(model@S), rownames(model@v0)))
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

##' Determine if the trajectory is sparse.
##' @noRd
is_trajectory_sparse <- function(x)
{
    if (identical(dim(x), c(0L, 0L)))
        return(FALSE)
    TRUE
}

##' Extract data in the internal matrix format
##'
##' @param m simulated data to extract.
##' @param ac available compartments in the simulated data.
##' @param sc selected compartments to extract from the simulated data
##'     and include in the matrix.
##' @param i subset of nodes to extract data from. If NULL, all
##'     available nodes are included.
##' @noRd
trajectory_as_is <- function(m, ac, sc, i)
{
    if (is.null(i)) {
        if (length(sc) == length(ac))
            return(m)
        i <- seq_len(nrow(m) %/% length(ac))
    }

    ## Extract subset of data.
    sc <- sort(match(sc, ac))
    i <- rep(sc, length(i)) + rep((i - 1) * length(ac), each = length(sc))
    m[i, , drop = FALSE]
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
trajectory <- function(model, compartments = NULL, node = NULL, as.is = FALSE)
{
    check_model_argument(model)

    if (is_trajectory_empty(model)) {
        stop("Please run the model first, the trajectory is empty.",
             call. = FALSE)
    }

    if (is(compartments, "formula"))
        compartments <- parse_formula(model, compartments)

    compartments <- match_compartments(model, compartments, as.is)

    ## Check the 'node' argument
    node <- check_node_argument(model, node)

    ## Check to extract data in internal matrix format
    if (isTRUE(as.is)) {
        if (!is.null(compartments$V)) {
            if (is_trajectory_sparse(model@V_sparse))
                return(trajectory_as_is(model@V_sparse, rownames(model@v0),
                                        compartments$V, node))

            return(trajectory_as_is(model@V, rownames(model@v0),
                                    compartments$V, node))
        }

        if (is_trajectory_sparse(model@U_sparse))
            return(trajectory_as_is(model@U_sparse, rownames(model@S),
                                    compartments$U, node))

        return(trajectory_as_is(model@U, rownames(model@S),
                                compartments$U, node))
    }

    ## Coerce the dense/sparse 'U' and 'V' matrices to a data.frame
    ## with one row per node and time-point with data from the
    ## specified discrete and continuous states.
    if (is_trajectory_sparse(model@U_sparse)) {
        Um <- model@U_sparse
    } else {
        Um <- model@U
    }
    Um_i <- match(compartments$U, rownames(model@S))
    Um_lbl <- rownames(model@S)

    if (is_trajectory_sparse(model@V_sparse)) {
        Vm <- model@V_sparse
    } else {
        Vm <- model@V
    }
    Vm_i <- match(compartments$V, rownames(model@v0))
    Vm_lbl <- rownames(model@v0)

    .Call(SimInf_trajectory,
          Um, Um_i, Um_lbl,
          Vm, Vm_i, Vm_lbl,
          model@tspan, Nn(model), node)
}
