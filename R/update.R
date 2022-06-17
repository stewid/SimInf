## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2022 Stefan Widgren
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

##' Update the initial compartment state u0 in each node
##'
##' @param model The model to update the initial compartment state
##'     \code{u0}.
##' @param value A \code{data.frame} with the initial state in each
##'     node. Each row is one node, and the number of rows in
##'     \code{u0} must match the number of nodes in \code{model}. Only
##'     the columns in \code{u0} with a name that matches a
##'     compartment in the \code{model} will be used.
##' @export
##' @examples
##' ## Create an SIR model object.
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the SIR model and plot the result.
##' set.seed(22)
##' result <- run(model)
##' plot(result)
##'
##' ## Update u0 and run the model again
##' update_u0(model) <- data.frame(S = 990, I = 10, R = 0)
##' result <- run(model)
##' plot(result)
setGeneric(
    "update_u0<-",
    signature = "model",
    function(model, value) {
        standardGeneric("update_u0<-")
    }
)

##' @rdname update_u0-set
##' @export
setMethod(
    "update_u0<-",
    signature(model = "SimInf_model"),
    function(model, value) {
        compartments <- rownames(model@S)
        value <- check_u0(value, compartments)
        if (!identical(nrow(value), n_nodes(model))) {
            stop("The number of rows in 'u0' must match nodes in 'model'.",
                 call. = FALSE)
        }
        model@u0 <- init_x0(value)
        model
    }
)

##' Update the initial continuous state v0 in each node
##'
##' @param model The model to update the initial continuous state
##'     \code{v0}.
##' @param v0 A \code{data.frame} with the initial continuosu state in
##'     each node. Each row is one node, and the number of rows in
##'     \code{v0} must match the number of nodes in \code{model}. Only
##'     the columns in \code{v0} with a name that matches a continuous
##'     state in \code{v0} in the \code{model} will be used.
##' @return a \code{SimInf_model} with the updated initial continuous
##'     state \code{v0}.
##' @export
##' @examples
##' ## Create an 'SISe' model with no infected individuals and no
##' ## infectious pressure (phi = 0, epsilon = 0).
##' model <- SISe(u0 = data.frame(S = 100, I = 0), tspan = 1:100,
##'               phi = 0, upsilon = 0.02, gamma = 0.1, alpha = 1,
##'               epsilon = 0, beta_t1 = 0.15, beta_t2 = 0.15,
##'               beta_t3 = 0.15, beta_t4 = 0.15, end_t1 = 91,
##'               end_t2 = 182, end_t3 = 273, end_t4 = 365)
##'
##' ## Run the 'SISe' model and plot the result.
##' set.seed(22)
##' result <- run(model)
##' plot(result)
##'
##' ## Update the infectious pressure 'phi' in 'v0' and run
##' ## the model again.
##' model <- update_v0(model, data.frame(phi = 1))
##' result <- run(model)
##' plot(result)
setGeneric(
    "update_v0",
    signature = "model",
    function(model, v0) {
        standardGeneric("update_v0")
    }
)

##' @rdname update_v0
##' @export
setMethod(
    "update_v0",
    signature(model = "SimInf_model"),
    function(model, v0) {
        variables <- rownames(model@v0)
        if (is.null(variables))
            variables <- character(0)

        v0 <- check_v0(v0, variables)
        if (!identical(nrow(v0), n_nodes(model))) {
            stop("The number of rows in 'v0' must match nodes in 'model'.",
                 call. = FALSE)
        }

        if (!identical(dim(model@v0), c(0L, 0L)))
            model@v0 <- init_x0(v0, "double")

        model
    }
)
