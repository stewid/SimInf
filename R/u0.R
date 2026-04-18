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

u0_target <- function(u0, target) {
    if (is.null(target))
        return(u0)

    if (target %in% c("SISe3", "SISe3_sp")) {
        u0$I_1 <- 0L
        u0$I_2 <- 0L
        u0$I_3 <- 0L
        return(u0)
    }

    if (target %in% c("SIS", "SISe", "SISe_sp")) {
        colnames(u0) <- c("key", "node", "S")
        u0$I <- 0L
        return(u0)
    }

    if (target %in% c("SIR")) {
        colnames(u0) <- c("key", "node", "S")
        u0$I <- 0L
        u0$R <- 0L
        return(u0)
    }

    if (target %in% c("SEIR")) {
        colnames(u0) <- c("key", "node", "S")
        u0$E <- 0L
        u0$I <- 0L
        u0$R <- 0L
        return(u0)
    }

    stop("Invalid 'target' for 'u0'.", call. = FALSE)
}

##' @rdname u0
##' @param time Only used when object is of class
##'     \code{SimInf_individual_events} object. The time-point that
##'     will be used to create u0. If left empty (the default), the
##'     earliest time among the events will be used.
##' @param target Only used when object is of class
##'     \code{SimInf_individual_events} object. The SimInf model
##'     ('SEIR', 'SIR', 'SIS', 'SISe3', 'SISe3_sp', 'SISe', or
##'     'SISe_sp') to target the events and u0 for. The default,
##'     \code{NULL}, creates an \code{u0}, but where the compartments
##'     might have to be renamed and post-processed to fit the
##'     specific use case.
##' @param age Only used when object is of class
##'     \code{SimInf_individual_events} object. An integer vector with
##'     break points in days for the ageing events. The default,
##'     \code{NULL}, creates an \code{u0} where all individuals belong
##'     to the same age category.
##' @export
setMethod(
    "u0",
    signature(object = "SimInf_individual_events"),
    function(object, time = NULL, target = NULL, age = NULL) {
        age <- check_age(age)
        target <- check_target(target, age)

        ## Determine the location and age for all individuals.
        individuals <- get_individuals(object, time)

        ## Ensure all nodes are included in u0.
        all_nodes <- unique(c(object@node, object@dest))
        all_nodes <- all_nodes[!is.na(all_nodes)]
        missing_nodes <- setdiff(all_nodes, individuals$node)
        S_columns <- paste0("S_", seq_along(age))

        if (nrow(individuals)) {
            ## Determine the age categories.
            age_category <- paste0("S_", findInterval(individuals$age, age))
            age_category <- c(age_category,
                              rep(NA_character_, length(missing_nodes)))

            ## Create u0.
            nodes <- c(individuals$node, missing_nodes)
            u0 <- as.data.frame.matrix(table(nodes, age_category))

            ## Ensure all age categories exist in u0
            age_category <- setdiff(S_columns, colnames(u0))
            if (length(age_category)) {
                u0 <- cbind(u0,
                            matrix(data = 0L,
                                   nrow = length(all_nodes),
                                   ncol = length(age_category),
                                   dimnames = list(NULL, age_category)))
            }
        } else {
            ## Create an empty u0.
            u0 <- as.data.frame.matrix(
                matrix(data = 0L,
                       nrow = length(all_nodes),
                       ncol = length(age),
                       dimnames = list(all_nodes, S_columns)))
        }

        u0 <- cbind(key = rownames(u0), u0)
        mode(u0$key) <- mode(all_nodes)
        rownames(u0) <- NULL
        u0 <- u0[order(u0$key), ]
        u0$node <- seq_len(nrow(u0))
        u0 <- u0[, c("key", "node", S_columns), drop = FALSE]

        u0_target(u0, target)
    }
)

##' Update the initial compartment state (\code{u0}) in each node
##'
##' Replace the initial state vector (\code{u0}) of a
##' \code{SimInf_model} object with a new \code{data.frame}. This
##' allows you to modify the starting conditions of a model without
##' recreating the object.
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
