## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2020 Stefan Widgren
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

##' Determine the compartments in the formula item and split
##' 'compartment1 + compartment2 + ...'. Moreover, trim whitespace.
##' @noRd
parse_formula_item <- function(x) {
    x <- unique(trimws(unlist(strsplit(x, "+", fixed = TRUE))))
    if (!length(x))
        stop("No compartments in formula specification.", call. = FALSE)
    x
}

define_dot <- function(a, b, args, ok_combine) {
    dot <- NULL

    if (!is.null(b)) {
        dot <- unlist(lapply(args, function(x) {
            if (any(b %in% x))
                return(x)
            NULL
        }), use.names = FALSE)
    }

    if (is.null(dot)) {
        if (isTRUE(ok_combine)) {
            dot <- unlist(args, use.names = FALSE)
        } else {
            dot <- args[[1]]
        }
    }

    dot
}

##' Replace '.' with compartments in the model.
##' @noRd
replace_dot <- function(a, b, args, ok_combine) {
    if (is.null(a))
        return(NULL)

    dot <- define_dot(a, b, args, ok_combine)
    a <- unlist(sapply(a, function(x) {
        if (identical(x, "."))
            x <- dot
        x
    }))

    unique(a)
}

parse_formula <- function(compartments, args, ok_combine) {
    compartments <- as.character(compartments)
    if (identical(length(compartments), 2L)) {
        lhs <- NULL
        rhs <- compartments[2]
    } else if (identical(length(compartments), 3L)) {
        lhs <- parse_formula_item(compartments[2])
        rhs <- compartments[3]
    } else {
        stop("Invalid formula specification of 'compartments'.", call. = FALSE)
    }

    ## Check if the rhs of the formula contains a condition.
    condition <- NULL
    if (any(regexpr("|", rhs, fixed = TRUE) > 1)) {
        condition <- sub("(^[^|]+)([|]?)(.*)$", "\\3", rhs)
        condition <- trimws(condition[nchar(condition) > 0])
        if (length(condition) != 1) {
            stop("Invalid formula specification of 'condition'.",
                 call. = FALSE)
        }

        rhs <- sub("(^[^|]+)([|]?)(.*)$", "\\1", rhs)
    }

    rhs <- parse_formula_item(rhs)
    lhs <- replace_dot(lhs, rhs, args, ok_combine)
    rhs <- replace_dot(rhs, lhs, args, ok_combine)

    list(lhs = lhs, rhs = rhs, condition = condition)
}

select_compartments <- function(x, from) {
    if (is.null(x))
        return(NULL)

    lapply(from, function(y) {
        x[x %in% y]
    })
}

##' Convert the compartment names to a named vector of indices.
##' Additionally, store all available compartments as an attribute.
##' @noRd
transform_compartments <- function(compartments, args) {
    if (is.null(compartments))
        return(NULL)

    mapply(function(a, b) {
        x <- match(a, b)
        names(x) <- b[x]
        attr(x, "available_compartments") <- b
        x
    }, compartments, args, SIMPLIFY = FALSE)
}

stop_if_combined_data <- function(lhs, rhs) {
    msg <- "Cannot combine data from different slots."

    if (!is.null(lhs)) {
        i <- vapply(lhs, function(x) {
            length(x) > 0
        }, logical(1))
        lhs <- names(lhs)[i]
        if (length(lhs) > 1)
            stop(msg, call. = FALSE)
    }

    if (!is.null(rhs)) {
        i <- vapply(rhs, function(x) {
            length(x) > 0
        }, logical(1))
        rhs <- names(rhs)[i]
        if (length(rhs) > 1)
            stop(msg, call. = FALSE)
    }

    if (!is.null(lhs) && !is.null(rhs) && lhs != rhs)
        stop(msg, call. = FALSE)

    invisible(NULL)
}

check_matched_data <- function(ok_combine, ok_lhs, lhs, rhs, compartments) {
    compartments <- setdiff(unlist(c(compartments$lhs, compartments$rhs)),
                            unlist(c(lhs, rhs)))
    if (length(compartments) > 0) {
        stop("Non-existing compartment(s) in model: ",
             paste0("'", compartments, "'", collapse = ", "),
             ".", call. = FALSE)
    }

    if (!isTRUE(ok_lhs) && !is.null(lhs))
        stop("Invalid formula specification of 'compartments'.", call. = FALSE)

    if (!isTRUE(ok_combine))
        stop_if_combined_data(lhs, rhs)

    invisible(NULL)
}

##' Match the selected 'compartments' argument in a function with the
##' available compartments in a model.
##'
##' @param compartments the names of the compartments to extract data
##'     from. The compartments can be specified as a formula or as a
##'     character vector.
##' @param ok_combine logical to indicate whether data from differnt
##'     slots can be combined.
##' @param ok_lhs logical to indicate whether a left-hand-side of a
##'     formula is allowed.
##' @param ... character vectors with available compartment names in
##'
##' the model.
##' @return a list with indices to the compartments in the available
##'     data-structures in the model.
##' @noRd
match_compartments <- function(compartments, ok_combine, ok_lhs, ...) {
    args <- list(...)

    if (methods::is(compartments, "formula")) {
        compartments <- parse_formula(compartments, args, ok_combine)
    } else {
        compartments <- list(lhs = NULL,
                             rhs = unique(trimws(as.character(compartments))),
                             condition = NULL)
    }

    lhs <- select_compartments(compartments$lhs, args)
    rhs <- select_compartments(compartments$rhs, args)
    condition <- compartments$condition

    check_matched_data(ok_combine, ok_lhs, lhs, rhs, compartments)

    ## If no compartments were selected, default to set 'rhs' to all
    ## available compartments.
    if (all(sapply(lhs, length) == 0, sapply(rhs, length) == 0))
        rhs <- args

    lhs <- transform_compartments(lhs, args)
    rhs <- transform_compartments(rhs, args)

    list(lhs = lhs, rhs = rhs, condition = condition)
}
