## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
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

##' Check integer arguments
##'
##' Raise an error if any of the arguments are non-integer.
##' @param len Expected length of the infectious pressure vector
##' @param ... The arguments to check
##' @return invisible(NULL)
##' @noRd
check_infectious_pressure_arg <- function(len, ...) {
    arg <- list(...)
    for (i in seq_len(length(arg))) {
        if (!is.numeric(arg[[i]]) || !is.null(dim(arg[[i]])) ||
            !identical(length(arg[[i]]), len)) {
            stop(paste0("Invalid '",
                        match.call(expand.dots = FALSE)$"..."[i],
                        "': must be numeric vector with length 'nrow(u0)'."),
                 call. = FALSE)
        }

        if (any(arg[[i]] < 0)) {
            stop(paste0("Invalid '",
                        match.call(expand.dots = FALSE)$"..."[i],
                        "': must be numeric vector with non-negative values."),
                 call. = FALSE)
        }
    }

    invisible(NULL)
}

##' Check integer arguments
##'
##' Raise an error if any of the arguments are non-integer.
##' @param ... The arguments to check
##' @return invisible(NULL)
##' @noRd
check_integer_arg <- function(...) {
    arg <- list(...)
    for (i in seq_len(length(arg))) {
        if (is.null(arg[[i]])) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$"..."[i],
                        "' is missing."),
                 call. = FALSE)
        }

        if (!is.numeric(arg[[i]]) ||
            any(is.na(arg[[i]])) ||
            !all(is_wholenumber(arg[[i]]))) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$"..."[i],
                        "' must be integer."),
                 call. = FALSE)
        }
    }

    invisible(NULL)
}

##' Check arguments for 'gdata'
##'
##' Raise an error if any of the arguments are not ok.
##' @param ... The arguments to check
##' @return invisible(NULL)
##' @noRd
check_gdata_arg <- function(...) {
    arg <- list(...)
    for (i in seq_len(length(arg))) {
        if (!is.numeric(arg[[i]]) || !identical(length(arg[[i]]), 1L)) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$"..."[i],
                        "' must be numeric of length 1."),
                 call. = FALSE)
        }
    }

    invisible(NULL)
}

##' Check arguments for 'ldata'
##'
##' Raise an error if any of the arguments are not ok.
##' @param len Exprected length of each data vector in '...'.
##' @param ... The arguments to check
##' @return invisible(NULL)
##' @noRd
check_ldata_arg <- function(len, ...) {
    arg <- list(...)
    for (i in seq_len(length(arg))) {
        if (!is.numeric(arg[[i]]) ||
            !is.atomic(arg[[i]]) ||
            (!identical(length(arg[[i]]), 1L) &&
             !identical(length(arg[[i]]), len))) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$"..."[i],
                        "' must be numeric of length 1 or 'nrow(u0)'."),
                 call. = FALSE)
        }
    }

    invisible(NULL)
}

##' Check arguments for interval endpoints
##'
##' Raise an error if any of the arguments are not ok.
##' @param len Exprected length of each of the interval endpoint
##' vectors
##' @param ... The arguments to check
##' @return invisible(NULL)
##' @noRd
check_end_t_arg <- function(len, ...) {
    arg <- list(...)
    names(arg) <- match.call(expand.dots = FALSE)$"..."

    for (i in seq_len(length(arg))) {
        if (!identical(length(arg[[i]]), len)) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$"..."[i],
                        "' must be of length 1 or 'nrow(u0)'."),
                 call. = FALSE)
        }
    }

    ## Check interval endpoints
    if (!all(0 <= arg$end_t1))
        stop("'end_t1' must be greater than or equal to '0'.", call. = FALSE)
    if (!all(arg$end_t1 < arg$end_t2))
        stop("'end_t1' must be less than 'end_t2'.", call. = FALSE)
    if (!all((arg$end_t4 < arg$end_t1) | (arg$end_t2 < arg$end_t3))) {
        stop(paste0("'end_t2' must be less than 'end_t3' or ",
                    "'end_t3' less than 'end_t1'."),
             call. = FALSE)
    }
    if (!all(arg$end_t3 < 364))
        stop("'end_t3' must be less than '364'.", call. = FALSE)
    if (!all(0 <= arg$end_t4))
        stop("'end_t4' must be greater than or equal to '0'.", call. = FALSE)
    if (!all(arg$end_t4 <= 365))
        stop("'end_t4' must be less than or equal to '365'.", call. = FALSE)
    if (!all((arg$end_t4 < arg$end_t1) | (arg$end_t3 < arg$end_t4))) {
        stop("'end_t4' must be less than 'end_t1' or greater than 'end_t3'.",
             call. = FALSE)
    }

    invisible(NULL)
}

##' Check model argument
##'
##' Raise an error if the model argument is not ok.
##' @param model the model to check.
##' @return invisible(NULL)
##' @noRd
check_model_argument <- function(model) {
    if (missing(model))
        stop("Missing 'model' argument.", call. = FALSE)
    if (!is(model, "SimInf_model"))
        stop("'model' argument is not a 'SimInf_model'.", call. = FALSE)

    invisible(NULL)
}

##' Check if wholenumbers
##'
##' Check that all values are wholenumbers, see example in integer {base}
##' @param x Value to check
##' @param tol Tolerance of the check
##' @return logical vector
##' @noRd
is_wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) {
    abs(x - round(x)) < tol
}

##' Check the node index argument
##'
##' Raise an error if the node argument is not ok.
##' @param model the model with nodes.
##' @param index the node index vector to check.
##' @return the node vector with unique nodes sorted in order.
##' @noRd
check_node_index_argument <- function(model, index) {
    if (is.null(index))
        return(NULL)

    if (!is.numeric(index) ||
        !all(is_wholenumber(index)) ||
        min(index) < 1 ||
        max(index) > n_nodes(model))
        stop("The node index must be an integer > 0 and <= number of nodes.",
             call. = FALSE)

    as.integer(sort(unique(index)))
}

##' Check the shift matrix 'N'
##'
##' Raise an error if the 'N' argument is not ok.
##' @param model the model with nodes.
##' @param N the shift matrix to check
##' @return the shift matrix.
##' @noRd
check_N <- function(N) {
    if (is.null(N))
        return(matrix(integer(0), nrow = 0, ncol = 0))

    if (!all(is.matrix(N), is.numeric(N)))
        stop("'N' must be an integer matrix.", call. = FALSE)

    if (!is.integer(N)) {
        if (!all(is_wholenumber(N)))
            stop("'N' must be an integer matrix.", call. = FALSE)
        storage.mode(N) <- "integer"
    }

    N
}

##' Check compartments
##'
##' Raise an error if the 'compartments' arguments is invalid.
##' @param compartments a character vector with the compartment names
##'     to check.
##' @return invisible(NULL)
##' @noRd
check_compartments <- function(compartments) {
    if (!is.atomic(compartments) || !is.character(compartments) ||
        !identical(compartments, make.names(compartments, unique = TRUE))) {
        stop("'compartments' must be specified in a character vector.",
             call. = FALSE)
    }

    invisible(NULL)
}

##' Check u0
##'
##' Raise an error if any of the 'u0' or 'compartments' arguments are
##' invalid.
##' @param u0 the initial state for the model.
##' @param compartments the compartments in u0.
##' @return u0 with columns ordered by the compartments.
##' @noRd
check_u0 <- function(u0, compartments) {
    check_compartments(compartments)

    ## Check u0
    if (!is.data.frame(u0))
        u0 <- as.data.frame(u0)
    if (!all(compartments %in% names(u0)))
        stop("Missing columns in u0.", call. = FALSE)

    u0[, compartments, drop = FALSE]
}

##' Check v0
##'
##' Raise an error if any of the 'v0' or 'compartments' arguments are
##' invalid.
##' @param v0 the initial continuous state for the model.
##' @param variables the variables in v0.
##' @return v0 with columns ordered by the variables.
##' @noRd
check_v0 <- function(v0, variables) {
    check_compartments(variables)

    ## Check v0
    if (!is.data.frame(v0))
        v0 <- as.data.frame(v0)
    if (!all(variables %in% names(v0)))
        stop("Missing columns in 'v0'.", call. = FALSE)

    v0[, variables, drop = FALSE]
}

##' Check distance matrix
##'
##' Raise an error if the distance argument is not ok.
##' @param distance The distance matrix between neighboring nodes
##' @return invisible(NULL)
##' @noRd
check_distance_matrix <- function(distance) {
    if (is.null(distance))
        stop("'distance' is missing.", call. = FALSE)
    if (!is(distance, "dgCMatrix")) {
        stop("The 'distance' argument must be of type 'dgCMatrix'.",
             call. = FALSE)
    }
    if (any(distance < 0))
        stop("All values in the 'distance' matrix must be >= 0.", call. = FALSE)

    invisible(NULL)
}

##' Check a package name
##'
##' From the "Writing R Extensions" manual: The mandatory ‘Package’
##' field gives the name of the package. This should contain only
##' (ASCII) letters, numbers and dot, have at least two characters and
##' start with a letter and not end in a dot.
##' @param name Character string with the package name.
##' @return invisible(NULL)
##' @noRd
check_package_name <- function(name) {
    pattern <- paste0("^", .standard_regexps()$valid_package_name, "$")

    if (any(is.null(name), !is.character(name), length(name) != 1,
            nchar(name) == 0, !grepl(pattern, name))) {
        stop("Malformed package name.", call. = FALSE)
    }

    invisible(NULL)
}

##' Check raw events identifiers
##'
##' Check that the identifiers used for 'id', 'node', and 'dest' are
##' valid, else raise an error.
##' @param ... The arguments ('id', 'node', and 'dest') to check.
##' @return invisible(NULL)
##' @noRd
check_raw_events_identifier <- function(...) {
    arg <- list(...)
    names(arg) <- match.call(expand.dots = FALSE)$"..."

    for (i in seq_len(length(arg))) {
        is_ok <- FALSE

        if (is.numeric(arg[[i]])) {
            is_ok <- all(is_wholenumber(arg[[i]]))
        } else {
            is_ok <- is.character(arg[[i]])
        }

        if (isTRUE(is_ok)) {
            is_ok <- !any(is.na(arg[[i]]))
        }

        if (!isTRUE(is_ok)) {
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$"..."[i],
                        "' must be an integer or character vector ",
                        "with non-NA values."),
                 call. = FALSE)
        }
    }

    invisible(NULL)
}

##' Check raw events
##'
##' Raise an error if any of the columns in the events data.frame are
##' invalid.
##' @param events a data.frame with raw events.
##' @return a raw events data.frame with the columns 'id', 'event',
##'     'time', 'node', and 'dest'.
##' @noRd
check_raw_events <- function(events) {
    columns <- c("id", "event", "time", "node", "dest")
    if (!is.data.frame(events))
        events <- as.data.frame(events)
    if (!all(columns %in% names(events)))
        stop("Missing columns in 'events'.", call. = FALSE)

    check_raw_events_identifier(events$id, events$node, events$dest)

    ## Do we have to recode the event type as a numerical value
    if (any(is.character(events$event), is.factor(events$event))) {
        event_names <- c("enter", "exit", "extTrans", "intTrans")
        if (!all(events$event %in% event_names)) {
            stop(paste0("'event' type must be 'enter', 'exit', ",
                        "'extTrans' or 'intTrans'."),
                 call. = FALSE)
        }

        ## Find indices to 'enter', 'internal transfer' and 'external
        ## transfer' events.
        i_enter <- which(events$event == "enter")
        i_intTrans <- which(events$event == "intTrans")
        i_extTrans <- which(events$event == "extTrans")

        ## Replace the character event type with a numerical value.
        events$.event <- rep(0L, nrow(events))
        events$.event[i_enter] <- 1L
        events$.event[i_intTrans] <- 2L
        events$.event[i_extTrans] <- 3L
    }

    events[, columns, drop = FALSE]
}
