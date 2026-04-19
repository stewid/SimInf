## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
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

##' Get the initial compartment state (\code{u0}) in each node
##'
##' Extract the initial state vector (\code{u0}) from a
##' \code{SimInf_model} object as a \code{data.frame}.
##'
##' The returned data frame has one row per node and one column per
##' compartment (e.g., \code{S}, \code{I}, \code{R}). This format is
##' convenient for inspection, modification, or exporting the initial
##' conditions.
##'
##' @param object A \code{SimInf_model} object.
##' @param ... Additional arguments.
##' @return a \code{data.frame} with the initial compartment state.
##' @export
##' @examples
##' ## Create an SIR model object.
##' model <- SIR(
##'   u0 = data.frame(S = 99, I = 1, R = 0),
##'   tspan = 1:100,
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Get the initial compartment state.
##' u0(model)
##'
##' ## Modify the initial state (e.g., add 10 infected individuals to
##' ## node 1).
##' new_u0 <- u0(model)
##' new_u0[1, "I"] <- new_u0[1, "I"] + 10
##'
##' ## Create a new model with the modified initial state.
##' new_model <- SIR(
##'   u0 = new_u0,
##'   tspan = 1:100,
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Alternatively, update the existing model using the setter:
##' u0(model) <- new_u0
## nolint start: brace_linter
setGeneric(
    "u0",
    signature = "object",
    function(object,
             ...)
        standardGeneric("u0")
)
## nolint end

##' @rdname u0
##' @export
setMethod(
    "u0",
    signature(object = "SimInf_model"),
    function(object, ...) {
        as.data.frame(t(object@u0))
    }
)

##' Update the initial compartment state (\code{u0}) in each node
##'
##' Replace the initial state vector (\code{u0}) of a
##' \code{SimInf_model} object with new data. This allows you to
##' modify the starting conditions of a model without recreating the
##' object.
##'
##' The \code{value} argument accepts a \code{data.frame},
##' \code{matrix}, or \code{named numeric vector}. If the input is not
##' a \code{data.frame}, it will be automatically coerced to one. The
##' function handles the following formats:
##' \itemize{
##'   \item \strong{Single Node}: If \code{value} is a named vector or
##'     a one-row matrix/data.frame, it is applied to the single node
##'     in the model.
##'   \item \strong{Multiple Nodes}: If \code{value} is a matrix or
##'     data.frame with multiple rows, each row corresponds to one
##'     node. The number of rows must exactly match the number of
##'     nodes in the \code{model}.
##'   \item \strong{Column Matching}: Column names must match the
##'     compartment names defined in the model (e.g., \code{"S"},
##'     \code{"I"}, \code{"R"}).  Only matching columns are used;
##'     extra columns are ignored, and missing compartments will
##'     trigger an error.
##' }
##'
##' The function validates the input and ensures the new state is
##' consistent with the model structure before updating.
##'
##' @param model A \code{SimInf_model} object.
##' @param value An object containing the new initial state. Can be a
##'     \code{data.frame}, \code{matrix}, or \code{named numeric
##'     vector}.  Non-data.frame inputs will be coerced to a
##'     \code{data.frame}.
##' @return The modified \code{SimInf_model} object.
##' @seealso \code{\link{v0<-}} for updating the initial continuous
##'     state.
##' @export
##' @examples
##' ## For reproducibility, set the seed.
##' set.seed(22)
##'
##' ## Create a single-node SIR model.
##' model <- SIR(
##'   u0 = data.frame(
##'     S = 99,
##'     I = 1,
##'     R = 0
##'   ),
##'   tspan = 1:100,
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Update u0 using a named vector (automatically coerced to one
##' ## row).
##' u0(model) <- c(
##'   S = 990,
##'   I = 10,
##'   R = 0
##' )
##'
##' result <- run(model)
##' plot(result)
##'
##' ## Create a multi-node model (2 nodes).
##' model_multi <- SIR(
##'   u0 = data.frame(
##'     S = c(100, 50),
##'     I = c(1, 0),
##'     R = c(0, 0)
##'   ),
##'   tspan = 1:100,
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Update u0 using a data.frame with multiple rows.
##' u0(model_multi) <- data.frame(
##'   S = c(200, 100),
##'   I = c(5, 2),
##'   R = c(0, 0)
##' )
##'
##' result <- run(model_multi)
##' plot(result)
## nolint start: brace_linter
setGeneric(
    "u0<-",
    signature = "model",
    function(model,
             value)
        standardGeneric("u0<-")
)
## nolint end

##' @rdname u0-set
##' @export
setMethod(
    "u0<-",
    signature(model = "SimInf_model"),
    function(model, value) {
        compartments <- rownames(model@S)
        value <- check_initial_state(value, compartments)
        if (!identical(nrow(value), n_nodes(model))) {
            stop("The number of rows in 'u0' must match nodes in 'model'.",
                 call. = FALSE)
        }
        model@u0 <- init_x0(value)
        methods::validObject(model)
        model
    }
)
