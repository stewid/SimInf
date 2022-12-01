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
##' u0(model) <- data.frame(S = 990, I = 10, R = 0)
##' result <- run(model)
##' plot(result)
setGeneric(
    "u0<-",
    signature = "model",
    function(model, value) {
        standardGeneric("u0<-")
    }
)

##' @rdname u0-set
##' @export
setMethod(
    "u0<-",
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
