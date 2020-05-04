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

## Determine the compartments in the formula item and split
## 'compartment1 + compartment2 + ...'. Moreover, trim whitespace and
## replace '.' with all compartments in the model.
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
    i <- !(x %in% compartments)
    if (any(i)) {
        stop("Non-existing compartment(s) in model: ",
             paste0("'", x[i], "'", collapse = ", "),
             ".", call. = FALSE)
    }
    x
}

parse_formula <- function(x, compartments) {
    x <- as.character(x)
    if (!identical(length(x), 2L))
        stop("Invalid formula specification of 'compartments'.", call. = FALSE)

    parse_formula_item(x[2], compartments)
}
