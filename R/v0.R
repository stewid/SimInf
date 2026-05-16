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

##' Update the initial continuous state (\code{v0}) in each node
##'
##' Replace the initial continuous state vector (\code{v0}) of a
##' \code{SimInf_model} object with new data. This allows you to
##' modify the starting conditions of continuous variables (e.g.,
##' environmental pathogen concentration) without recreating the
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
##'     continuous state variable names defined in the model (e.g.,
##'     \code{"phi"}).  Only matching columns are used; extra columns
##'     are ignored, and missing variables will trigger an error.
##' }
##'
##' The function validates the input and ensures the new state is
##' consistent with the model structure before updating.
##'
##' @param model A \code{SimInf_model} object.
##' @param value An object containing the new initial continuous
##'     state.  Can be a \code{data.frame}, \code{matrix}, or
##'     \code{named numeric vector}.  Non-data.frame inputs will be
##'     coerced to a \code{data.frame}.
##' @return The modified \code{SimInf_model} object.
##' @seealso \code{\link{u0<-}} for updating the initial discrete
##'     compartment state.
##' @export
##' @examples
##' ## For reproducibility, set the seed.
##' set.seed(22)
##'
##' ## Create an 'SISe' model with no infected individuals and no
##' ## infectious pressure (phi = 0).
##' model <- SISe(
##'   u0 = data.frame(S = 100, I = 0),
##'   tspan = 1:100,
##'   phi = 0,
##'   upsilon = 0.02,
##'   gamma = 0.1,
##'   alpha = 1,
##'   epsilon = 0,
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
##' ## Run the 'SISe' model and plot the result.
##' result <- run(model)
##' plot(result)
##'
##' ## Update the infectious pressure 'phi' in 'v0' using a named
##' ## vector.  (Automatically coerced to one row for the single
##' ## node).
##' v0(model) <- c(phi = 1)
##' result <- run(model)
##' plot(result)
##'
##' ## For a multi-node model, use a data.frame with multiple rows:
##' ## v0(model) <- data.frame(phi = c(1.0, 0.5, 0.0))
## nolint start: brace_linter
setGeneric(
    "v0<-",
    signature = "model",
    function(model,
             value)
        standardGeneric("v0<-")
)
## nolint end

##' @rdname v0-set
##' @export
setMethod(
    "v0<-",
    signature(model = "SimInf_model"),
    function(model, value) {
        variables <- rownames(model@v0)
        if (is.null(variables))
            variables <- character(0)

        value <- check_initial_state(value, variables)
        if (!identical(nrow(value), n_nodes(model))) {
            stop("The number of rows in 'v0' must match nodes in 'model'.",
                 call. = FALSE)
        }

        if (!identical(dim(model@v0), c(0L, 0L)))
            model@v0 <- init_x0(value, "double")

        methods::validObject(model)
        model
    }
)
