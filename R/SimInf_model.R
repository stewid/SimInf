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

##' Check if a SimInf_model object is valid
##'
##' @param object The SimInf_model_object to check.
##' @include classes.R
##' @include valid.R
##' @noRd
valid_SimInf_model_object <- function(object) {
    ## Check events
    methods::validObject(object@events)

    errors <- c(valid_tspan(object),
                valid_replicates(object),
                valid_u0(object),
                valid_U(object),
                valid_v0(object),
                valid_V(object),
                valid_S(object),
                valid_G(object),
                valid_ldata(object),
                valid_gdata(object))

    if (length(errors))
        return(errors)
    TRUE
}

## Assign the function as the validity method for the class.
setValidity("SimInf_model", valid_SimInf_model_object)

##' Create a \code{SimInf_model} object
##'
##' Construct a low-level \code{SimInf_model} object. This function is
##' typically used internally by model constructors (e.g.,
##' \code{SIR()}, \code{mparse()}) or for advanced usage where custom
##' model definitions (e.g., user-provided C code or non-standard
##' matrices) are required.
##'
##' @param G \strong{Dependency Graph}.  Indicates which transition
##'     rates need updating after a state transition.  Can be provided
##'     as a sparse matrix (class \code{dgCMatrix}) or a dense matrix.
##'     If a dense matrix is provided, it is automatically converted
##'     to a sparse format internally.  See
##'     \code{\linkS4class{SimInf_model}} for detailed matrix layout.
##'
##' @param S \strong{State Transition Matrix}.  Defines the change in
##'     the state vector for each transition.  Can be provided as a
##'     sparse matrix (class \code{dgCMatrix}) or a dense matrix.  If
##'     a dense matrix is provided, it is automatically converted to a
##'     sparse format internally.  See
##'     \code{\linkS4class{SimInf_model}} for detailed matrix layout.
##'
##' @param U \strong{Result Matrix} (integer matrix).  Usually empty
##'     at creation.  See \code{\linkS4class{SimInf_model}} for
##'     detailed matrix layout.
##'
##' @param ldata \strong{Local Data}.
##'     Parameters specific to each node. Can be:
##'     \itemize{
##'       \item A \code{data.frame} with one row per node.
##'       \item A matrix where each column \code{ldata[, j]} is the
##'       data vector for node \code{j}.
##'     }
##'     Passed to transition rate and post-step functions.
##'
##' @param gdata \strong{Global Data} (numeric vector).  Parameters
##'     common to all nodes. Passed to transition rate and post-step
##'     functions.
##'
##' @param tspan \strong{Time Span} (numeric or Date vector).
##'     Increasing time points for output. If \code{Date}, converted
##'     to days with names, where \code{tspan[1]} becomes the day of
##'     the year of the first year of \code{tspan}. The dates are
##'     added as names to the numeric vector.
##'
##' @param u0 \strong{Initial State}.  Initial number of individuals
##'     per compartment/node. Can be:
##'     \itemize{
##'       \item A matrix (\eqn{N_c \times N_n}).
##'       \item A \code{data.frame} with columns corresponding to
##'       compartments.
##'       \item Any object coercible to a \code{data.frame} (e.g., a
##'         named numeric vector will be coerced to a one-row
##'         \code{data.frame}).
##'     }
##'
##' @param events \strong{Scheduled Events}.  A \code{data.frame}
##'     defining the event schedule (see
##'     \code{\linkS4class{SimInf_events}}).
##'
##' @param V \strong{Continuous State Result Matrix} (numeric matrix).
##'     Usually empty at creation.  See
##'     \code{\linkS4class{SimInf_model}} for layout.
##'
##' @param v0 \strong{Initial Continuous State} (numeric matrix).
##'     Initial values for continuous states per node.
##'
##' @param E \strong{Select Matrix} (matrix or \code{data.frame}).
##'     Defines which compartments are affected by events and their
##'     sampling weights.
##'     \itemize{
##'       \item \strong{Matrix}: Standard sparse matrix.
##'       \item \strong{data.frame}: Must have columns
##'         \code{compartment} and \code{select}.  Optional column
##'         \code{value} (default \code{1}) sets the weight.
##'     }
##'     See \code{\linkS4class{SimInf_events}} for usage details.
##'
##' @param N \strong{Shift Matrix} (matrix or \code{data.frame}).
##'     Defines how individuals are moved between compartments during
##'     events.
##'     \itemize{
##'       \item \strong{Matrix}: Standard integer matrix.
##'       \item \strong{data.frame}: Must have columns
##'         \code{compartment}, \code{shift}, and \code{value}
##'         (integer offset).
##'     }
##'     See \code{\linkS4class{SimInf_events}} for usage details.
##'
##' @param C_code \strong{C Source Code} (character vector).  Optional
##'     C code for custom transition rates. If provided, it is
##'     compiled and loaded when \code{run()} is called.
##'
##' @return A \code{\linkS4class{SimInf_model}} object.
##'
##' @seealso \code{\link{SIR}}, \code{\link{SEIR}}, \code{\link{SIS}},
##'     \code{\link{SISe}} for examples of compartment model
##'     constructors that handle argument validation and matrix setup.
##'     \code{\link{mparse}} for creating custom models using a simple
##'     string syntax.  \code{\linkS4class{SimInf_model}} for details
##'     on the class structure and slots.  \code{\link{SimInf_events}}
##'     for details on the event schedule format.
##'
##' @include init.R
##' @export
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
                         C_code = NULL) {
    u0 <- init_x0(u0)
    G <- init_sparse_matrix(G)
    S <- init_sparse_matrix(S)
    ldata <- init_data_matrix(ldata)
    gdata <- init_data_vector(gdata)
    U <- init_output_matrix(U)
    v0 <- init_x0(v0, "double", TRUE)
    V <- init_output_matrix(V, "double")
    C_code <- init_C_code(C_code)
    tspan <- init_tspan(tspan)

    ## Check events
    if (!any(is.null(events), is.data.frame(events)))
        stop("'events' must be NULL or a data.frame.", call. = FALSE)
    if (is.data.frame(E))
        E <- E_from_data_frame(E, rownames(S))
    if (is.data.frame(N))
        N <- N_from_data_frame(N, rownames(S))
    events <- SimInf_events(E = E, N = N, events = events, t0 = tspan$t0)

    methods::new("SimInf_model",
                 G          = G,
                 S          = S,
                 U          = U,
                 ldata      = ldata,
                 gdata      = gdata,
                 tspan      = tspan$tspan,
                 u0         = u0,
                 v0         = v0,
                 V          = V,
                 events     = events,
                 replicates = 1L,
                 C_code     = C_code)
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
## nolint start: brace_linter
setGeneric(
    "gdata",
    signature = "model",
    function(model)
        standardGeneric("gdata")
)
## nolint end

##' @rdname gdata
##' @export
setMethod(
    "gdata",
    signature(model = "SimInf_model"),
    function(model) {
        model@gdata
    }
)

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
## nolint start: brace_linter
setGeneric(
    "gdata<-",
    signature = "model",
    function(model,
             parameter,
             value)
        standardGeneric("gdata<-")
)
## nolint end

##' @rdname gdata-set
##' @export
setMethod(
    "gdata<-",
    signature(model = "SimInf_model"),
    function(model, parameter, value) {
        ## Check parameter argument
        if (missing(parameter))
            stop("Missing 'parameter' argument.", call. = FALSE)
        if (!is.character(parameter))
            stop("'parameter' argument must be a character.", call. = FALSE)

        ## Check value argument
        if (missing(value))
            stop("Missing 'value' argument.", call. = FALSE)
        if (!is.numeric(value))
            stop("'value' argument must be a numeric.", call. = FALSE)

        model@gdata[parameter] <- value

        methods::validObject(model)
        model
    }
)

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
## nolint start: brace_linter
setGeneric(
    "ldata",
    signature = "model",
    function(model,
             node)
        standardGeneric("ldata")
)
## nolint end

##' @rdname ldata
##' @export
setMethod(
    "ldata",
    signature(model = "SimInf_model"),
    function(model, node) {
        ## Check node argument
        if (missing(node))
            stop("Missing 'node' argument.", call. = FALSE)
        if (!is.numeric(node) || !identical(length(node), 1L) || node < 1)
            stop("Invalid 'node' argument.", call. = FALSE)

        model@ldata[, node]
    }
)
