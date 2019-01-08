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

##' Class \code{"SimInf_model"}
##'
##' Class to handle the siminf data model
##' @section Slots:
##' \describe{
##'   \item{G}{
##'     Dependency graph that indicates the transition rates that need
##'     to be updated after a given state transition has occured.
##'     A non-zero entry in element \code{G[i, i]} indicates that transition
##'     rate \code{i} needs to be recalculated if the state transition
##'     \code{j} occurs. Sparse matrix (\eqn{Nt \times Nt}) of object class
##'     \code{\linkS4class{dgCMatrix}}.
##'   }
##'   \item{S}{
##'     Each column corresponds to a state transition, and execution
##'     of state transition \code{j} amounts to adding the \code{S[,
##'     j]} column to the state vector \code{u[, i]} of node \emph{i}
##'     where the transition occurred. Sparse matrix (\eqn{Nc \times
##'     Nt}) of object class \code{\linkS4class{dgCMatrix}}.
##'   }
##'   \item{U}{
##'     The result matrix with the number of individuals in each
##'     compartment in every node. \code{U[, j]} contains the number
##'     of individuals in each compartment at
##'     \code{tspan[j]}. \code{U[1:Nc, j]} contains the number of
##'     individuals in node 1 at \code{tspan[j]}. \code{U[(Nc + 1):(2
##'     * Nc), j]} contains the number of individuals in node 2 at
##'     \code{tspan[j]} etc. Integer matrix (\eqn{N_n N_c \times}
##'     \code{length(tspan)}).
##'   }
##'   \item{U_sparse}{
##'     If the model was run to write the solution to a sparse matrix
##'     (\code{dgCMatrix}) the \code{U_sparse} contains the data and
##'     \code{U} is empty. The layout of the data in \code{U_sparse}
##'     is identical to \code{U}. Please note that \code{U_sparse}
##'     is numeric and \code{U} is integer.
##'   }
##'   \item{V}{
##'     The result matrix for the real-valued continuous
##'     state. \code{V[, j]} contains the real-valued state of the
##'     system at \code{tspan[j]}. Numeric matrix
##'     (\eqn{N_n}\code{dim(ldata)[1]} \eqn{\times}
##'     \code{length(tspan)}).
##'   }
##'   \item{V_sparse}{
##'     If the model was run to write the solution to a sparse matrix
##'     (\code{dgCMatrix}) the \code{V_sparse} contains the data and
##'     \code{V} is empty. The layout of the data in \code{V_sparse}
##'     is identical to \code{V}.
##'   }
##'   \item{ldata}{
##'     A matrix with local data for the nodes. The column \code{ldata[, j]}
##'     contains the local data vector for the node \code{j}. The local
##'     data vector is passed as an argument to the transition rate
##'     functions and the post time step function.
##'   }
##'   \item{gdata}{
##'     A numeric vector with global data that is common to all nodes.
##'     The global data vector is passed as an argument to the
##'     transition rate functions and the post time step function.
##'   }
##'   \item{tspan}{
##'     A vector of increasing time points where the state of each node is
##'     to be returned.
##'   }
##'   \item{u0}{
##'     The initial state vector (\eqn{N_c \times N_n}) with
##'     the number of individuals in each compartment in every node.
##'   }
##'   \item{v0}{
##'      The initial value for the real-valued continuous state.
##'      Numeric matrix (\code{dim(ldata)[1]} \eqn{\times N_n}).
##'   }
##'   \item{events}{
##'     Scheduled events \code{\linkS4class{SimInf_events}}
##'   }
##'   \item{C_code}{
##'     Character vector with optional model C code. If non-empty, the
##'     C code is written to a temporary C-file when the \code{run}
##'     method is called.  The temporary C-file is compiled and the
##'     resulting DLL is dynamically loaded. The DLL is unloaded and
##'     the temporary files are removed after running the model.
##'   }
##' }
##' @include SimInf_events.R
##' @export
##' @importFrom methods validObject
##' @importClassesFrom Matrix dgCMatrix
setClass("SimInf_model",
         slots = c(G        = "dgCMatrix",
                   S        = "dgCMatrix",
                   U        = "matrix",
                   U_sparse = "dgCMatrix",
                   ldata    = "matrix",
                   gdata    = "numeric",
                   tspan    = "numeric",
                   u0       = "matrix",
                   V        = "matrix",
                   V_sparse = "dgCMatrix",
                   v0       = "matrix",
                   events   = "SimInf_events",
                   C_code   = "character"),
         validity = function(object) {
             ## Check events
             errors <- validObject(object@events)
             if (!isTRUE(errors))
                 return(errors)

             ## Check tspan.
             if (!is.double(object@tspan)) {
                 return("Input time-span must be a double vector.")
             } else if (any(length(object@tspan) < 2,
                            any(diff(object@tspan) <= 0),
                            any(is.na(object@tspan)))) {
                 return("Input time-span must be an increasing vector.")
             }

             ## Check u0.
             if (!identical(storage.mode(object@u0), "integer"))
                 return("Initial state 'u0' must be an integer matrix.")
             if (any(object@u0 < 0L))
                 return("Initial state 'u0' has negative elements.")
             Nn_u0 <- dim(object@u0)[2]

             ## Check U.
             if (!identical(storage.mode(object@U), "integer"))
                 return("Output state 'U' must be an integer matrix.")
             if (any(object@U < 0L))
                 return("Output state 'U' has negative elements.")

             ## Check v0.
             if (!identical(storage.mode(object@v0), "double"))
                 return("Initial model state 'v0' must be a double matrix.")
             if ((dim(object@v0)[1] > 0)) {
                 r <- rownames(object@v0)
                 if (is.null(r) || any(nchar(r) == 0))
                     return("'v0' must have rownames")
                 if (!identical(dim(object@v0)[2], Nn_u0))
                     return("The number of nodes in 'u0' and 'v0' must match.")
             }

             ## Check V.
             if (!identical(storage.mode(object@V), "double"))
                 return("Output model state 'V' must be a double matrix.")

             ## Check S.
             if (!all(is_wholenumber(object@S@x)))
                 return("'S' matrix must be an integer matrix.")

             ## Check that S and events@E have identical compartments
             if ((dim(object@S)[1] > 0) && (dim(object@events@E)[1] > 0)) {
                 if (!identical(rownames(object@S), rownames(object@events@E)))
                     return("'S' and 'E' must have identical compartments")
             }

             ## Check G.
             Nt <- dim(object@S)[2]
             if (!identical(dim(object@G), c(Nt, Nt)))
                 return("Wrong size of dependency graph.")

             ## Check that transitions exist in G.
             transitions <- rownames(object@G)
             if (is.null(transitions))
                 return("'G' must have rownames that specify transitions.")
             transitions <- sub("^[[:space:]]*", "", sub("[[:space:]]*$", "", transitions))
             if (!all(nchar(transitions) > 0))
                 return("'G' must have rownames that specify transitions.")

             ## Check that the format of transitions are valid.
             ## "X1 + X2 + ... + Xn -> Y1 + Y2 + ... + Yn"
             ## is expected, where X2, ..., Xn and Y2, ..., Yn are optional.
             transitions <- strsplit(transitions, split = "->", fixed = TRUE)
             if (!all(sapply(transitions, length) == 2))
                 return("'G' rownames have invalid transitions.")

             ## Check that transitions and S have identical compartments
             transitions <- unlist(transitions)
             transitions <- unlist(strsplit(transitions, split = "+", fixed = TRUE))
             transitions <- sub("^[[:space:]]*", "", sub("[[:space:]]*$", "", transitions))
             transitions <- unique(transitions)
             transitions <- transitions[transitions != "@"]
             transitions <- sub("^[[:digit:]]+[*]", "", transitions)
             if (!all(transitions %in% rownames(object@S)))
                 return("'G' and 'S' must have identical compartments")

             ## Check ldata.
             if (!is.double(object@ldata))
                 return("'ldata' matrix must be a double matrix.")
             Nn_ldata <- dim(object@ldata)[2]
             if (Nn_ldata > 0 && !identical(Nn_ldata, Nn_u0))
                 return("The number of nodes in 'u0' and 'ldata' must match.")

             ## Check gdata.
             if (!is.double(object@gdata))
                 return("'gdata' must be a double vector.")

             TRUE
         }
)

## Utility function to coerce the data.frame to a transposed matrix.
as_t_matrix <- function(x) {
    n_col <- ncol(x)
    n_row <- nrow(x)
    lbl <- colnames(x)
    x <- t(data.matrix(x))
    attributes(x) <- NULL
    dim(x) <- c(n_col, n_row)
    rownames(x) <- lbl
    x
}

##' Create a \code{SimInf_model}
##'
##' @template G-param
##' @template S-param
##' @template U-param
##' @template ldata-param
##' @template gdata-param
##' @template tspan-param
##' @param u0 The initial state vector. Either a matrix (\eqn{N_c
##'     \times N_n}) or a a \code{data.frame} with the number of
##'     individuals in each compartment in every node.
##' @param events A \code{data.frame} with the scheduled events.
##' @param V The result matrix for the real-valued continous
##'     compartment state (\eqn{N_n}\code{dim(ldata)[1]} \eqn{\times}
##'     \code{length(tspan)}).  \code{V[, j]} contains the real-valued
##'     state of the system at \code{tspan[j]}.
##' @param v0 The initial continuous state vector in every node.
##'     (\code{dim(ldata)[1]} \eqn{\times N_N}). The continuous state
##'     vector is updated by the specific model during the simulation
##'     in the post time step function.
##' @param E Sparse matrix to handle scheduled events, see
##'     \code{\linkS4class{SimInf_events}}.
##' @param N Sparse matrix to handle scheduled events, see
##'     \code{\linkS4class{SimInf_events}}.
##' @param C_code Character vector with optional model C code. If
##'     non-empty, the C code is written to a temporary C-file when
##'     the \code{run} method is called.  The temporary C-file is
##'     compiled and the resulting DLL is dynamically loaded. The DLL
##'     is unloaded and the temporary files are removed after running
##'     the model.
##' @return \linkS4class{SimInf_model}
##' @export
##' @importFrom methods as
##' @importFrom methods is
##' @importFrom methods new
SimInf_model <- function(G,
                         S,
                         tspan,
                         events = NULL,
                         ldata  = NULL,
                         gdata  = NULL,
                         U      = NULL,
                         u0     = NULL,
                         v0     = NULL,
                         V      = NULL,
                         E      = NULL,
                         N      = NULL,
                         C_code = NULL)
{
    ## Check u0
    if (is.null(u0))
        stop("'u0' is NULL")
    if (is.data.frame(u0))
        u0 <- as_t_matrix(u0)
    if (!all(is.matrix(u0), is.numeric(u0)))
        stop("u0 must be an integer matrix")
    if (!is.integer(u0)) {
        if (!all(is_wholenumber(u0)))
            stop("u0 must be an integer matrix")
        storage.mode(u0) <- "integer"
    }

    ## Check G
    if (!is.null(G)) {
        if (!is(G, "dgCMatrix"))
            G <- as(G, "dgCMatrix")
    }

    ## Check S
    if (!is.null(S)) {
        if (!is(S, "dgCMatrix"))
            S <- as(S, "dgCMatrix")
    }

    ## Check ldata
    if (is.null(ldata))
        ldata <- matrix(numeric(0), nrow = 0, ncol = 0)
    if (is.data.frame(ldata))
        ldata <- as_t_matrix(ldata)

    ## Check gdata
    if (is.null(gdata))
        gdata <- numeric(0)
    if (is.data.frame(gdata)) {
        if (!identical(nrow(gdata), 1L))
            stop("When 'gdata' is a data.frame, it must have one row.")
        gdata <- unlist(gdata)
    }

    ## Check U
    if (is.null(U)) {
        U <- matrix(integer(0), nrow = 0, ncol = 0)
    } else {
        if (!is.integer(U)) {
            if (!all(is_wholenumber(U)))
                stop("U must be an integer")
            storage.mode(U) <- "integer"
        }

        if (!is.matrix(U)) {
            if (!identical(length(U), 0L))
                stop("U must be equal to 0 x 0 matrix")
            dim(U) <- c(0, 0)
        }
    }

    ## Check v0
    if (is.null(v0)) {
        v0 <- matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
        if (is.data.frame(v0))
            v0 <- as_t_matrix(v0)
        if (!all(is.matrix(v0), is.numeric(v0)))
            stop("v0 must be a numeric matrix")

        if (!identical(storage.mode(v0), "double"))
            storage.mode(v0) <- "double"
    }

    ## Check V
    if (is.null(V)) {
        V <- matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
        if (!is.numeric(V))
            stop("V must be numeric")

        if (!identical(storage.mode(V), "double"))
            storage.mode(V) <- "double"

        if (!is.matrix(V)) {
            if (!identical(length(V), 0L))
                stop("V must be equal to 0 x 0 matrix")
            dim(V) <- c(0, 0)
        }
    }

    ## Check tspan
    if (is(tspan, "Date")) {
        ## Coerce the date vector to a numeric vector as days, where
        ## tspan[1] becomes the day of the year of the first year of
        ## the tspan date vector. The dates are added as names to the
        ## numeric vector.
        t0 <- as.numeric(as.Date(format(tspan[1], "%Y-01-01"))) - 1
        tspan_lbl <- format(tspan, "%Y-%m-%d")
        tspan <- as.numeric(tspan) - t0
        names(tspan) <- tspan_lbl
    } else {
        t0 <- NULL
    }
    storage.mode(tspan) <- "double"

    ## Check events
    if (!any(is.null(events), is.data.frame(events)))
        stop("'events' must be NULL or a data.frame")
    events <- SimInf_events(E = E, N = N, events = events, t0 = t0)

    ## Check C code
    if (is.null(C_code))
        C_code <- character(0)

    new("SimInf_model",
        G      = G,
        S      = S,
        U      = U,
        ldata  = ldata,
        gdata  = gdata,
        tspan  = tspan,
        u0     = u0,
        v0     = v0,
        V      = V,
        events = events,
        C_code = C_code)
}

##' Set a template for where to write the U result matrix
##'
##' Using a sparse U result matrix can save a lot of memory if the
##' model contains many nodes and time-points, but where only a few of
##' the data points are of interest for post-processing.
##'
##' Using a sparse U result matrix can save a lot of memory if the
##' model contains many nodes and time-points, but where only a few of
##' the data points are of interest for post-processing. To use this
##' feature, a template has to be defined for which data points to
##' record. This is done using a \code{data.frame} that specifies the
##' time-points (column \sQuote{time}) and nodes (column
##' \sQuote{node}) to record the state of the compartments, see
##' \sQuote{Examples}. The specified time-points, nodes and
##' compartments must exist in the model, or an error is raised. Note
##' that specifying a template only affects which data-points are
##' recorded for post-processing, it does not affect how the solver
##' simulates the trajectory.
##' @param model The \code{model} to set a template for the result
##'     matrix \code{U}.
##' @param value A \code{data.frame} that specify the nodes,
##'     time-points and compartments to record the number of
##'     individuals at \code{tspan}. Use \code{NULL} to reset the
##'     model to record the number of inidividuals in each compartment
##'     in every node at each time-point in tspan.
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
##' U(model) <- df
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
##' U(model) <- df
##' result <- run(model, threads = 1)
##' trajectory(result)
##'
##' ## It is possible to use an empty 'data.frame' to specify
##' ## that no data-points should be recorded for the trajectory.
##' U(model) <- data.frame()
##' result <- run(model, threads = 1)
##' trajectory(result)
##'
##' ## Use 'NULL' to reset the model to record data for every node at
##' ## each time-point in tspan.
##' U(model) <- NULL
##' result <- run(model, threads = 1)
##' trajectory(result)
"U<-" <- function(model, value)
{
    check_model_argument(model)

    if (!is.null(value)) {
        if (!is.data.frame(value))
            stop("'value' argument is not a 'data.frame'")

        if (nrow(value) > 0) {
            ## Sort the data.frame by time and node.
            value <- value[order(value$time, value$node),
                           c("time", "node", rownames(model@S))]

            ## Match nodes and for each matched node create an index
            ## to all of its compartments in the U matrix.
            i <- match(value$node, seq_len(Nn(model)))
            if (any(is.na(i)))
                stop("Unable to match all nodes")
            i <- rep((i - 1) * Nc(model), each = Nc(model)) + seq_len(Nc(model))

            ## Match time-points to tspan and repeat each time-point
            ## for every compartment in the model.
            j <- match(value$time, model@tspan)
            if (any(is.na(j)))
                stop("Unable to match all time-points to tspan")
            j <- rep(j, each = Nc(model))

            ## Coerce the compartments part of the data.frame to a
            ## logical vector that match the rows of compartments in
            ## the U matrix.
            value <- as.logical(t(as.matrix(value[, -(1:2)])))
            value[is.na(value)] <- FALSE

            ## Keep only compartments and time-points that are marked
            ## with TRUE.
            i <- i[value]
            j <- j[value]
        } else {
            i <- numeric(0)
            j <- numeric(0)
        }

        ## Specify dimension.
        d <- c(Nn(model) * Nc(model), length(model@tspan))

        ## Clear dense result matrix.
        model@U = matrix(integer(0), nrow = 0, ncol = 0)
    } else {
        ## Clear sparse result matrix.
        i <- numeric(0)
        j <- numeric(0)
        d <- c(0, 0)
    }

    ## Create sparse template.
    model@U_sparse <- as(sparseMatrix(i, j, dims = d), "dgCMatrix")

    model
}

##' Set a template for where to write the V result matrix
##'
##' Using a sparse V result matrix can save a lot of memory if the
##' model contains many nodes and time-points, but where only a few of
##' the data points are of interest for post-processing.
##'
##' Using a sparse V result matrix can save a lot of memory if the
##' model contains many nodes and time-points, but where only a few of
##' the data points are of interest for post-processing. To use this
##' feature, a template has to be defined for which data points to
##' record. This is done using a \code{data.frame} that specifies the
##' time-points (column \sQuote{time}) and nodes (column
##' \sQuote{node}) to record the state of the continuous state
##' compartments, see \sQuote{Examples}. The specified time-points,
##' nodes and compartments must exist in the model, or an error is
##' raised. Note that specifying a template only affects which
##' data-points are recorded for post-processing, it does not affect
##' how the solver simulates the trajectory.
##' @param model The \code{model} to set a template for the result
##'     matrix \code{V}.
##' @param value A \code{data.frame} that specify the nodes,
##'     time-points and compartments of when to record the real-valued
##'     continuous state at \code{tspan}. Use \code{NULL} to reset the
##'     model to record the real-valued continuous state in every node
##'     at each time-point in tspan.
##' @export
##' @importFrom methods as
##' @importFrom Matrix sparseMatrix
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
##' result <- run(model, threads = 1)
##'
##' ## Display the continuous state variable 'phi' for every node at
##' ## each time-point in tspan.
##' trajectory(result, compartments = "phi")
##'
##' ## Assume we are only interested in nodes '2' and '4' at the
##' ## time-points '3' and '5'
##' df <- data.frame(time = c(3, 5, 3, 5),
##'                  node = c(2, 2, 4, 4),
##'                  phi = c(TRUE, TRUE, TRUE, TRUE))
##' V(model) <- df
##' result <- run(model, threads = 1)
##' trajectory(result, compartments = "phi")
##'
##' ## It is possible to use an empty 'data.frame' to specify
##' ## that no data-points should be recorded for the trajectory.
##' V(model) <- data.frame()
##' result <- run(model, threads = 1)
##' trajectory(result, compartments = "phi")
##'
##' ## Use 'NULL' to reset the model to record data for every node at
##' ## each time-point in tspan.
##' V(model) <- NULL
##' result <- run(model, threads = 1)
##' trajectory(result, compartments = "phi")
"V<-" <- function(model, value)
{
    check_model_argument(model)

    if (!is.null(value)) {
        if (!is.data.frame(value))
            stop("'value' argument is not a 'data.frame'")

        if (nrow(value) > 0) {
            ## Sort the data.frame by time and node.
            value <- value[order(value$time, value$node),
                           c("time", "node", rownames(model@v0))]

            ## Match nodes and for each matched node create an index
            ## to all of its continuous state compartments in the V
            ## matrix.
            i <- match(value$node, seq_len(Nn(model)))
            if (any(is.na(i)))
                stop("Unable to match all nodes")
            i <- rep((i - 1) * Nd(model), each = Nd(model)) + seq_len(Nd(model))

            ## Match time-points to tspan and repeat each time-point
            ## for every continuous state compartment in the model.
            j <- match(value$time, model@tspan)
            if (any(is.na(j)))
                stop("Unable to match all time-points to tspan")
            j <- rep(j, each = Nd(model))

            ## Coerce the compartments part of the data.frame to a
            ## logical vector that match the rows of continuous state
            ## compartments in the V matrix.
            value <- as.logical(t(as.matrix(value[, -(1:2)])))
            value[is.na(value)] <- FALSE

            ## Keep only the compartments and time-points that are
            ## marked with TRUE
            i <- i[value]
            j <- j[value]
        } else {
            i <- numeric(0)
            j <- numeric(0)
        }

        ## Specify dimension.
        d <- c(Nn(model) * Nd(model), length(model@tspan))

        ## Clear dense result matrix
        model@V <- matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
        ## Clear sparse result matrix
        i <- numeric(0)
        j <- numeric(0)
        d <- c(0, 0)
    }

    ## Create sparse template
    model@V_sparse <- as(sparseMatrix(i = i, j = j, dims = d), "dgCMatrix")

    model
}

##' Extract number of nodes in a model
##'
##' Extract number of nodes in a model.
##' @param model the \code{model} object to extract the number of
##'     nodes from.
##' @return the number of nodes in the model.
##' @export
##' @examples
##' ## Create an 'SIR' model with 100 nodes, with 99 susceptible,
##' ## 1 infected and 0 recovered in each node.
##' u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Display the number of nodes in the model.
##' Nn(model)
Nn <- function(model)
{
    check_model_argument(model)
    dim(model@u0)[2]
}

## Number of compartments
Nc <- function(model)
{
    check_model_argument(model)
    dim(model@S)[1]
}

## Number of transitions
Nt <- function(model)
{
    check_model_argument(model)
    dim(model@G)[1]
}

## Number of continuous state variables
Nd <- function(model)
{
    check_model_argument(model)
    dim(model@v0)[1]
}

## Use 'R CMD SHLIB' to compile the C code for the model.
do_compile_model <- function(filename)
{
    ## Include directive for "SimInf.h"
    include <- system.file("include", package = "SimInf")
    Sys.setenv(PKG_CPPFLAGS=sprintf("-I%s", shQuote(include)))

    ## Compile the model C code using the running version of R.
    wd <- setwd(dirname(filename))
    cmd <- paste(shQuote(file.path(R.home(component="bin"), "R")),
                 "CMD SHLIB",
                 shQuote(paste0(basename(filename), ".c")))
    compiled <- system(cmd, intern = TRUE)
    setwd(wd)

    lib <- paste0(filename, .Platform$dynlib.ext)
    if (!file.exists(lib))
        stop(compiled)

    lib
}

## Check if model contains C code
contains_C_code <- function(model)
{
    if (nchar(paste0(model@C_code, collapse = "\n")))
        return(TRUE)
    FALSE
}

##' Run the SimInf stochastic simulation algorithm
##'
##' @param model The siminf model to run.
##' @param threads Number of threads. Default is NULL, i.e. to use all
##'     available processors.
##' @param solver Which numerical solver to utilize. Default is 'ssm'.
##' @return \code{\link{SimInf_model}} object with result from simulation.
##' @references \itemize{
##'   \item Bauer P, Engblom S, Widgren S
##'   (2016) "Fast Event-Based Epidemiological Simulations on National Scales"
##'   International Journal of High Performance Computing
##'   Applications, 30(4), 438-453. doi:10.1177/1094342016635723
##'
##'   \item Bauer P., Engblom S. (2015) Sensitivity Estimation and
##'   Inverse Problems in Spatial Stochastic Models of Chemical
##'   Kinetics. In: Abdulle A., Deparis S., Kressner D., Nobile F.,
##'   Picasso M. (eds) Numerical Mathematics and Advanced Applications
##'   - ENUMATH 2013. Lecture Notes in Computational Science and
##'   Engineering, vol 103. Springer, Cham. Doi:
##'   10.1007/978-3-319-10705-9_51
##' }
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
##' result <- run(model, threads = 1)
##'
##' ## Plot the proportion of susceptible, infected and recovered
##' ## individuals.
##' plot(result)
setGeneric("run",
           signature = "model",
           function(model,
                    threads = NULL,
                    solver  = c("ssm", "aem"))
               standardGeneric("run"))

##' @rdname run
##' @export
##' @importFrom methods validObject
setMethod("run",
          signature(model = "SimInf_model"),
          function(model, threads, solver)
          {
              solver <- match.arg(solver)

              ## Check that SimInf_model contains all data structures
              ## required by the siminf solver and that they make sense
              validObject(model);

              if (contains_C_code(model)) {
                  ## Write the C code to a temporary file
                  filename <- tempfile("SimInf-")
                  on.exit(unlink(paste0(filename,
                                        c(".c", ".o", .Platform$dynlib.ex))))
                  writeLines(model@C_code, con = paste0(filename, ".c"))

                  lib <- do_compile_model(filename)
                  dll <- dyn.load(lib)
                  on.exit(dyn.unload(lib), add = TRUE)

                  ## Create expression to parse
                  expr <- ".Call(dll$SimInf_model_run, model, threads, solver)"
              } else {
                  ## The model name
                  name <- as.character(class(model))

                  ## The model C run function
                  run_fn <- paste0(name, "_run")

                  ## Create expression to parse
                  expr <- ".Call(run_fn, model, threads, solver, PACKAGE = 'SimInf')"
              }

              ## Run the model. Re-throw any error without the call
              ## included in the error message to make it cleaner.
              tryCatch(eval(parse(text = expr)), error = function(e) {
                  stop(e$message, call. = FALSE)
              })
          }
)

##' @importFrom stats quantile
##' @noRd
summary_U <- function(object)
{
    cat("Compartments\n")
    cat("------------\n")

    d <- dim(object@U)
    if (identical(d, c(0L, 0L)))
        d <- dim(object@U_sparse)
    if (identical(d, c(0L, 0L))) {
        cat(" - Empty, please run the model first\n")
    } else if (is.null(rownames(object@S))) {
        stop("Undefined model compartments")
    } else {
        qq <- lapply(rownames(object@S), function(compartment) {
            x <- as.numeric(trajectory(object, compartment, as.is = TRUE))
            qq <- quantile(x)
            qq <- c(qq[1L:3L], mean(x), qq[4L:5L])
        })
        qq <- do.call("rbind", qq)
        colnames(qq) <- c("Min.", "1st Qu.", "Median",
                          "Mean", "3rd Qu.", "Max.")
        rownames(qq) <- paste0(" ", rownames(object@S))
        print.table(qq, digits = 3)
    }
}

##' @importFrom stats quantile
##' @noRd
summary_V <- function(object)
{
    cat("Continuous state variables\n")
    cat("--------------------------\n")

    if (Nd(object) > 0) {
        d <- dim(object@V)
        if (identical(d, c(0L, 0L)))
            d <- dim(object@V_sparse)
        if (identical(d, c(0L, 0L))) {
            cat(" - Empty, please run the model first\n")
        } else if (is.null(rownames(object@v0))) {
            stop("Undefined continuous state variables")
        } else {
            qq <- lapply(rownames(object@v0), function(compartment) {
                x <- as.numeric(trajectory(object, compartment, as.is = TRUE))
                qq <- quantile(x)
                qq <- c(qq[1L:3L], mean(x), qq[4L:5L])
            })
            qq <- do.call("rbind", qq)
            colnames(qq) <- c("Min.", "1st Qu.", "Median",
                              "Mean", "3rd Qu.", "Max.")
            rownames(qq) <- rownames(object@v0)
            print.table(qq, digits = 3)
        }
    } else {
        cat(" - None\n")
    }
}

summary_ldata <- function(object)
{
    ## Local model parameters
    cat("Local data\n")
    cat("----------\n")

    if (dim(object@ldata)[1] > 0) {
        qq <- t(apply(object@ldata, 1, function(x) {
            qq <- quantile(x)
            c(qq[1L:3L], mean(x), qq[4L:5L])
        }))
        colnames(qq) <- c("Min.", "1st Qu.", "Median",
                          "Mean", "3rd Qu.", "Max.")
        rownames(qq) <- paste0(" ", rownames(object@ldata))
        print.table(qq, digits = 3)
    } else {
        cat(" - None\n")
    }
}

summary_gdata <- function(object)
{
    ## Global model parameters
    cat("Global data\n")
    cat("-----------\n")

    gdata <- data.frame(Parameter = names(object@gdata), Value = object@gdata)
    if (nrow(gdata) > 0) {
        print.data.frame(gdata, right = FALSE, row.names = FALSE)
    } else {
        cat(" - None\n")
    }
}

##' Determine in-degree for each node in a model
##'
##' The number of nodes with inward \emph{external transfer} events to
##' each node.
##' @param model determine in-degree for each node in the model.
##' @return vector with in-degree for each node.
##' @export
##' @examples
##' ## Create an 'SIR' model with 1600 nodes and initialize
##' ## it with example data.
##' model <- SIR(u0 = u0_SIR(), tspan = 1:1460, events = events_SIR(),
##'              beta   = 0.16, gamma  = 0.077)
##'
##' ## Display indegree for each node in the model.
##' plot(indegree(model))
indegree <- function(model)
{
    check_model_argument(model)

    ## Default indegree is 0
    id <- integer(Nn(model))

    ## Determine indegree from data
    i <- which(model@events@event == 3L)
    if (length(i) > 0) {
        idd <- tapply(model@events@node[i], model@events@dest[i],
                      function(x) {length(unique(x))})
        id[as.integer(dimnames(idd)[[1]])] <- idd
    }

    id
}

##' Determine out-degree for each node in a model
##'
##' The number nodes that are connected with \emph{external transfer}
##' events from each node.
##' @param model determine out-degree for each node in the model.
##' @return vector with out-degree for each node.
##' @export
##' @examples
##' ## Create an 'SIR' model with 1600 nodes and initialize
##' ## it with example data.
##' model <- SIR(u0 = u0_SIR(), tspan = 1:1460, events = events_SIR(),
##'              beta   = 0.16, gamma  = 0.077)
##'
##' ## Display outdegree for each node in the model.
##' plot(outdegree(model))
outdegree <- function(model)
{
    check_model_argument(model)

    ## Default outdegree is 0
    od <- integer(Nn(model))

    ## Determine oudegree from data
    i <- which(model@events@event == 3L)
    if (length(i) > 0) {
        odd <- tapply(model@events@dest[i], model@events@node[i],
                      function(x) {length(unique(x))})
        od[as.integer(dimnames(odd)[[1]])] <- odd
    }

    od
}

##' @importFrom stats quantile
##' @noRd
summary_events <- function(object)
{
    cat("Scheduled events\n")
    cat("----------------\n")

    if (length(object@events@event) > 0) {
        ## Summarise exit events
        i <- which(object@events@event == 0L)
        cat(sprintf(" Exit: %i\n", length(i)))

        ## Summarise enter events
        i <- which(object@events@event == 1L)
        cat(sprintf(" Enter: %i\n", length(i)))

        ## Summarise internal transfer events
        i <- which(object@events@event == 2L)
        cat(sprintf(" Internal transfer: %i\n", length(i)))

        ## Summarise external transfer events
        i <- which(object@events@event == 3L)
        cat(sprintf(" External transfer: %i\n", length(i)))

        if (length(i) > 0) {
            ## Summarise network
            cat("\nNetwork summary\n")
            cat("---------------\n")
            id <- indegree(object)
            od <- outdegree(object)
            qq_id <- quantile(id)
            qq_id <- c(qq_id[1L:3L], mean(id), qq_id[4L:5L])
            qq_od <- quantile(od)
            qq_od <- c(qq_od[1L:3L], mean(od), qq_od[4L:5L])
            qq <- rbind(qq_id, qq_od)
            colnames(qq) <- c("Min.", "1st Qu.", "Median",
                              "Mean", "3rd Qu.", "Max.")
            rownames(qq) <- c(" Indegree:", " Outdegree:")
            print.table(qq, digits = 3)
        }
    } else {
        cat(" - None\n")
    }
}

summary_transitions <- function(object)
{
    cat("Transitions\n")
    cat("-----------\n")

    cat(paste0(" ", rownames(object@G), collapse = "\n"), sep = "\n")
}

##' Brief summary of \code{SimInf_model}
##'
##' @param object The SimInf_model \code{object}
##' @return None (invisible 'NULL').
##' @export
##' @importFrom methods show
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
##' ## Brief summary of the model
##' model
##'
##' ## Run the model and save the result
##' result <- run(model, threads = 1)
##'
##' ## Brief summary of the result. Note that 'U' and 'V' are
##' ## non-empty after running the model.
##' result
setMethod("show",
          signature(object = "SimInf_model"),
          function (object)
          {
              ## The model name
              cat(sprintf("Model: %s\n", as.character(class(object))))
              cat(sprintf("Number of nodes: %i\n", Nn(object)))
              cat(sprintf("Number of transitions: %i\n", Nt(object)))
              show(object@events)

              cat("\n")
              summary_gdata(object)

              if (!is.null(rownames(object@ldata))) {
                  cat("\n")
                  summary_ldata(object)
              }

              if (Nd(object) > 0) {
                  cat("\n")
                  summary_V(object)
              }

              cat("\n")
              summary_U(object)

              invisible(object)
          }
)

##' Detailed summary of a \code{SimInf_model} object
##'
##' @param object The \code{SimInf_model} object
##' @param ... Additional arguments affecting the summary produced.
##' @return None (invisible 'NULL').
##' @export
setMethod("summary",
          signature(object = "SimInf_model"),
          function(object, ...)
          {
              ## The model name
              cat(sprintf("Model: %s\n", as.character(class(object))))

              ## Nodes
              cat(sprintf("Number of nodes: %i\n\n", Nn(object)))

              summary_transitions(object)

              cat("\n")
              summary_gdata(object)

              if (!is.null(rownames(object@ldata))) {
                  cat("\n")
                  summary_ldata(object)
              }

              cat("\n")
              summary_events(object)

              if (Nd(object) > 0) {
                  cat("\n")
                  summary_V(object)
              }

              cat("\n")
              summary_U(object)
          }
)

##' Extract the events from a \code{SimInf_model} object
##'
##' Extract the scheduled events from a \code{SimInf_model} object.
##' @param model The \code{model} to extract the events from.
##' @return \code{\linkS4class{SimInf_events}} object.
##' @export
##' @examples
##' ## Create an SIR model that includes scheduled events.
##' model <- SIR(u0     = u0_SIR(),
##'              tspan  = 1:(4 * 365),
##'              events = events_SIR(),
##'              beta   = 0.16,
##'              gamma  = 0.077)
##'
##' ## Extract the scheduled events from the model and display summary
##' summary(events(model))
##'
##' ## Extract the scheduled events from the model and plot them
##' plot(events(model))
events <- function(model)
{
    check_model_argument(model)
    model@events
}

##' Extract the shift matrix from a \code{SimInf_model} object
##'
##' Utility function to extract the shift matrix \code{events@@N} from
##' a \code{SimInf_model} object, see
##' \code{\linkS4class{SimInf_events}}
##' @param model The \code{model} to extract the shift matrix
##'     \code{events@@N} from.
##' @return A mtrix.
##' @export
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Extract the shift matrix from the model
##' shift_matrix(model)
shift_matrix <- function(model)
{
    check_model_argument(model)
    model@events@N
}

##' Set the shift matrix for a \code{SimInf_model} object
##'
##' Utility function to set \code{events@@N} in a \code{SimInf_model}
##' object, see \code{\linkS4class{SimInf_events}}
##' @param model The \code{model} to set the shift matrix
##'     \code{events@@N}.
##' @param value A matrix.
##' @return \code{SimInf_model} object
##' @export
##' @importFrom methods is
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Set the shift matrix
##' shift_matrix(model) <- matrix(c(2, 1, 0), nrow = 3)
##'
##' ## Extract the shift matrix from the model
##' shift_matrix(model)
"shift_matrix<-" <- function(model, value)
{
    check_model_argument(model)

    ## Check value
    if (is.null(value))
        value <- matrix(integer(0), nrow = 0, ncol = 0)
    if (!all(is.matrix(value), is.numeric(value)))
        stop("'value' must be an integer matrix")
    if (!is.integer(value)) {
        if (!all(is_wholenumber(value)))
            stop("'value' must be an integer matrix")
        storage.mode(value) <- "integer"
    }
    if (!identical(Nc(model), dim(value)[1]))
        stop("'value' must have one row for each compartment in the model")

    dimnames(value) <- list(rownames(model@events@E),
                            as.character(seq_len(dim(value)[2])))
    model@events@N <- value

    model
}

##' Extract the select matrix from a \code{SimInf_model} object
##'
##' Utility function to extract \code{events@@E} from a
##' \code{SimInf_model} object, see \code{\linkS4class{SimInf_events}}
##' @param model The \code{model} to extract the select matrix
##'     \code{E} from.
##' @return \code{\linkS4class{dgCMatrix}} object.
##' @export
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Extract the select matrix from the model
##' select_matrix(model)
select_matrix <- function(model)
{
    check_model_argument(model)
    model@events@E
}

##' Set the select matrix for a \code{SimInf_model} object
##'
##' Utility function to set \code{events@@E} in a \code{SimInf_model}
##' object, see \code{\linkS4class{SimInf_events}}
##' @param model The \code{model} to set the select matrix for.
##' @param value A matrix.
##' @export
##' @importFrom methods as
##' @importFrom methods is
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Set the select matrix
##' select_matrix(model) <- matrix(c(1, 0, 0, 1, 1, 1, 0, 0, 1), nrow = 3)
##'
##' ## Extract the select matrix from the model
##' select_matrix(model)
"select_matrix<-" <- function(model, value)
{
    check_model_argument(model)

    if (!is(value, "dgCMatrix"))
        value <- as(value, "dgCMatrix")

    if (!identical(Nc(model), dim(value)[1]))
        stop("'value' must have one row for each compartment in the model")

    dimnames(value) <- list(rownames(model@events@E),
                            as.character(seq_len(dim(value)[2])))
    model@events@E <- value

    model
}

##' Extract global data from a \code{SimInf_model} object
##'
##' The global data is a numeric vector that is common to all nodes.
##' The global data vector is passed as an argument to the transition
##' rate functions and the post time step function.
##' @param model The \code{model} to get global data from.
##' @return a numeric vector
##' @export
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Set 'beta' to a new value
##' gdata(model, "beta") <- 2
##'
##' ## Extract the global data vector that is common to all nodes
##' gdata(model)
gdata <- function(model)
{
    check_model_argument(model)
    model@gdata
}

##' Set a global data parameter for a \code{SimInf_model} object
##'
##' The global data is a numeric vector that is common to all nodes.
##' The global data vector is passed as an argument to the transition
##' rate functions and the post time step function.
##' @param model The \code{model} to set a global model parameter for.
##' @param parameter The name of the parameter to set.
##' @param value A numeric value.
##' @return a \code{SimInf_model} object
##' @export
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Set 'beta' to a new value
##' gdata(model, "beta") <- 2
##'
##' ## Extract the global data vector that is common to all nodes
##' gdata(model)
"gdata<-" <- function(model, parameter, value)
{
    check_model_argument(model)

    ## Check paramter argument
    if (missing(parameter))
        stop("Missing 'parameter' argument")
    if (!is.character(parameter))
        stop("'parameter' argument must be a character")

    ## Check value argument
    if (missing(value))
        stop("Missing 'value' argument")
    if (!is.numeric(value))
        stop("'value' argument must be a numeric")

    model@gdata[parameter] <- value

    model
}

##' Extract local data from a node
##'
##' The local data is a numeric vector that is specific to a node.
##' The local data vector is passed as an argument to the transition
##' rate functions and the post time step function.
##' @param model The \code{model} to get local data from.
##' @param node index to node to extract local data from.
##' @return a numeric vector
##' @export
##' @examples
##' ## Create an 'SISe' model with 1600 nodes.
##' model <- SISe(u0 = u0_SISe(), tspan = 1:100, events = events_SISe(),
##'               phi = 0, upsilon = 1.8e-2, gamma = 0.1, alpha = 1,
##'               beta_t1 = 1.0e-1, beta_t2 = 1.0e-1, beta_t3 = 1.25e-1,
##'               beta_t4 = 1.25e-1, end_t1 = c(91, 101), end_t2 = c(182, 185),
##'               end_t3 = c(273, 275), end_t4 = c(365, 360), epsilon = 0)
##'
##' ## Display local data from the first two nodes.
##' ldata(model, node = 1)
##' ldata(model, node = 2)
ldata <- function(model, node)
{
    check_model_argument(model)

    ## Check node argument
    if (missing(node))
        stop("Missing 'node' argument")
    if (!is.numeric(node) || !identical(length(node), 1L) || node < 1)
        stop("Invalid 'node' argument")

    model@ldata[, node]
}
