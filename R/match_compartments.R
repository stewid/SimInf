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

##' Match the 'compartments' argument in a function with the available
##' compartments in a model.
##'
##' @param compartments the names of the compartments to extract data
##'     from. The compartments can be specified as a formula or as a
##'     character vector.
##' @param ok_combine logical to indicate whether data from differnt
##'     slots can be combined.
##' @param ... character vectors with available compartment names in
##'     the model.
##' @return a list with indices to the compartments in the available
##'     data-structures in the model.
##' @noRd
match_compartments <- function(compartments, ok_combine, ...) {
    args <- list(...)

    if (is(compartments, "formula")) {
        compartments <- parse_formula(
            compartments, unlist(args, use.names = FALSE))
    }

    compartments <- unique(as.character(compartments))

    result <- lapply(args, function(x) {
        compartments[compartments %in% x]
    })

    compartments <- setdiff(compartments, unlist(result))
    if (length(compartments) > 0) {
        stop("Non-existing compartment(s) in model: ",
             paste0("'", compartments, "'", collapse = ", "),
             ".", call. = FALSE)
    }

    if (!isTRUE(ok_combine) && all(sapply(result, length))) {
        stop("Cannot combine data from different slots.",
             call. = FALSE)
    }

    if (all(sapply(result, length) == 0))
        result <- args

    mapply(match, result, args, SIMPLIFY = FALSE)
}
