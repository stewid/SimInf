## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2023 Stefan Widgren
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

##' Get the initial compartment state
##'
##' @param object The object to get the initial compartment state
##'     \code{u0} from.
##' @param ... Additional arguments.
##' @return a \code{data.frame} with the initial compartment state.
##' @export
##' @examples
##' ## Create an SIR model object.
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Get the initial compartment state.
##' u0(model)
setGeneric(
    "u0",
    signature = "object",
    function(object, ...) {
        standardGeneric("u0")
    }
)

##' @rdname u0
##' @export
setMethod(
    "u0",
    signature(object = "SimInf_model"),
    function(object, ...) {
        as.data.frame(t(object@u0))
    }
)

##' @rdname u0
##' @param time the time-point that will be used to create u0. If left
##'     empty (the default), the earliest time among the events will
##'     be used.
##' @param target the SimInf model ('SEIR', 'SIR', 'SIS', 'SISe3',
##'     'SISe3_sp', 'SISe', or 'SISe_sp') to target the events and u0
##'     for. The default, \code{NULL}, creates an \code{u0}, but where
##'     the compartments might have to be renamed and post-processed
##'     to fit the specific use case.
##' @param age an integer vector with break points in days for the
##'     ageing events. The default, \code{NULL}, creates an \code{u0}
##'     where all individuals belong to the same age category.
##' @export
setMethod(
    "u0",
    signature(object = "SimInf_indiv_events"),
    function(object, time = NULL, target = NULL, age = NULL) {
        ## Check for valid target model.
        if (!is.null(target)) {
            target <- match.arg(target, c("SEIR", "SIR", "SIS",
                                          "SISe3", "SISe3_sp", "SISe",
                                          "SISe_sp"))
        }

        if (is.null(age)) {
            age <- 0
        } else {
            ## Check for valid age.
            age <- sort(unique(as.integer(age)))
            if (any(age <= 0))
                stop("'age' must be an integer vector with values > 0.")
            age <- c(0, age)
        }

        ## Determine the location and age for all individuals.
        individuals <- get_individuals(object, time)

        ## Ensure all nodes are included in u0.
        all_nodes <- unique(c(object@node, object@dest))
        all_nodes <- all_nodes[!is.na(all_nodes)]
        all_nodes <- sort(all_nodes)
        missing_nodes <- setdiff(all_nodes, individuals$node)

        if (nrow(individuals)) {
            ## Determine the age categories.
            age_category <- paste0("S_", findInterval(individuals$age, age))
            age_category <- c(age_category,
                              rep(NA_character_, length(missing_nodes)))

            ## Create u0.
            nodes <- c(individuals$node, missing_nodes)
            u0 <- as.data.frame.matrix(table(nodes, age_category))

            ## Ensure all age categories exist in u0
            age_category <- setdiff(paste0("S_", seq_len(length(age))),
                                    colnames(u0))
            if (length(age_category)) {
                u0 <- cbind(u0,
                            matrix(data = 0L,
                                   nrow = length(all_nodes),
                                   ncol = length(age_category),
                                   dimnames = list(
                                       NULL,
                                       age_category)))
            }
        } else {
            ## Create an empty u0.
            u0 <- as.data.frame.matrix(
                matrix(data = 0L,
                       nrow = length(all_nodes),
                       ncol = length(age),
                       dimnames = list(
                           all_nodes,
                           paste0("S_", seq_len(length(age))))))
        }

        u0 <- u0[, paste0("S_", seq_len(length(age))), drop = FALSE]
        u0 <- cbind(key = rownames(u0),
                    node = seq_len(nrow(u0)),
                    u0)
        mode(u0$key) <- mode(all_nodes)
        rownames(u0) <- NULL

        u0
    }
)

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
        methods::validObject(model)
        model
    }
)
