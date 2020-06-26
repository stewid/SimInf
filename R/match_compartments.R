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
##' 'compartment1 + compartment2 + ...'. Moreover, trim whitespace and
##' replace '.' with all compartments in the model.
##' @noRd
parse_formula_item <- function(x, compartments) {
    x <- unlist(strsplit(x, "+", fixed = TRUE))
    x <- trimws(x)
    x <- unlist(sapply(x, function(y) {
        if (identical(y, "."))
            y <- compartments
        y
    }))
    x <- unique(as.character(x))
    if (!length(x))
        stop("No compartments in formula specification.", call. = FALSE)
    x
}

parse_formula <- function(x, compartments) {
    lhs <- NULL
    rhs <- NULL
    condition <- NULL

    x <- as.character(x)
    if (identical(length(x), 2L)) {
        rhs <- parse_formula_item(x[2], compartments)
    } else if (identical(length(x), 3L)) {
        lhs <- parse_formula_item(x[2], compartments)
        rhs <- parse_formula_item(x[3], compartments)

        ## Check if the rhs of the formula contains a condition.
        if (any(regexpr("|", rhs, fixed = TRUE) > 1)) {
            condition <- sub("(^[^|]+)([|]?)(.*)$", "\\3", rhs)
            condition <- trimws(condition[nchar(condition) > 0])
            if (length(condition) != 1) {
                stop("Invalid formula specification of 'condition'.",
                     call. = FALSE)
            }

            rhs <- sub("(^[^|]+)([|]?)(.*)$", "\\1", rhs)
            rhs <- parse_formula_item(rhs, compartments)
        }
    } else {
        stop("Invalid formula specification of 'compartments'.", call. = FALSE)
    }

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

    if (isTRUE(ok_combine))
        return(invisible(NULL))

    if (all(!is.null(lhs), (length(lhs) > 1), all(sapply(lhs, length))))
        stop("Cannot combine data from different slots.", call. = FALSE)
    if (all(!is.null(rhs), (length(rhs) > 1), all(sapply(rhs, length))))
        stop("Cannot combine data from different slots.", call. = FALSE)

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

    if (is(compartments, "formula")) {
        compartments <- parse_formula(
            compartments, unlist(args, use.names = FALSE))
    } else {
        compartments <- list(lhs = NULL,
                             rhs = unique(as.character(compartments)),
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
