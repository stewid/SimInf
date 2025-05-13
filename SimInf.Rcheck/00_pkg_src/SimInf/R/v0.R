## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2024 Stefan Widgren
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

##' Update the initial continuous state v0 in each node
##'
##' @param model The model to update the initial continuous state
##'     \code{v0}.
##' @param value the initial continuous state in each node. Must be a
##'     \code{data.frame} or an object that can be coerced to a
##'     \code{data.frame}. A named numeric vector will be coerced to a
##'     one-row \code{data.frame}. Each row is one node, and the
##'     number of rows in \code{v0} must match the number of nodes in
##'     \code{model}. Only the columns in \code{v0} with a name that
##'     matches a continuous state in \code{v0} in the \code{model}
##'     will be used
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
##' v0(model) <- data.frame(phi = 1)
##' result <- run(model)
##' plot(result)
setGeneric(
    "v0<-",
    signature = "model",
    function(model, value) {
        standardGeneric("v0<-")
    }
)

##' @rdname v0-set
##' @export
setMethod(
    "v0<-",
    signature(model = "SimInf_model"),
    function(model, value) {
        variables <- rownames(model@v0)
        if (is.null(variables))
            variables <- character(0)

        value <- check_v0(value, variables)
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
