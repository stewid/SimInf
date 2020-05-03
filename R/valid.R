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

valid_tspan <- function(object) {
    if (!is.double(object@tspan)) {
        return("Input time-span must be a double vector.")
    } else if (any(length(object@tspan) < 1,
                   any(diff(object@tspan) <= 0),
                   any(is.na(object@tspan)))) {
        return("Input time-span must be an increasing vector.")
    }

    character(0);
}

valid_u0 <- function(object) {
    if (!identical(storage.mode(object@u0), "integer"))
        return("Initial state 'u0' must be an integer matrix.")
    if (any(object@u0 < 0L))
        return("Initial state 'u0' has negative elements.")

    character(0);
}

valid_U <- function(object) {
    if (!identical(storage.mode(object@U), "integer"))
        return("Output state 'U' must be an integer matrix.")
    if (any(object@U < 0L) || any(object@U_sparse < 0, na.rm = TRUE))
        return("Output state 'U' has negative elements.")

    character(0);
}

valid_v0 <- function(object) {
    if (!identical(storage.mode(object@v0), "double"))
        return("Initial model state 'v0' must be a double matrix.")
    if ((dim(object@v0)[1] > 0)) {
        r <- rownames(object@v0)
        if (is.null(r) || any(nchar(r) == 0))
            return("'v0' must have rownames.")
        if (!identical(dim(object@v0)[2], dim(object@u0)[2]))
            return("The number of nodes in 'u0' and 'v0' must match.")
    }

    character(0);
}

valid_V <- function(object) {
    if (!identical(storage.mode(object@V), "double"))
        return("Output model state 'V' must be a double matrix.")

    character(0);
}

valid_S <- function(object) {
    if (!all(is_wholenumber(object@S@x)))
        return("'S' matrix must be an integer matrix.")

    ## Check that S and events@E have identical compartments
    if ((dim(object@S)[1] > 0) && (dim(object@events@E)[1] > 0)) {
        if (is.null(rownames(object@S)) || is.null(rownames(object@events@E)))
            return("'S' and 'E' must have rownames matching the compartments.")
        if (!identical(rownames(object@S), rownames(object@events@E)))
            return("'S' and 'E' must have identical compartments.")
    }

    character(0);
}

valid_G <- function(object) {
    Nt <- dim(object@S)[2]
    if (!identical(dim(object@G), c(Nt, Nt)))
        return("Wrong size of dependency graph.")

    ## Check that transitions exist in G.
    transitions <- rownames(object@G)
    if (is.null(transitions))
        return("'G' must have rownames that specify transitions.")
    transitions <- trimws(transitions)
    if (!all(nchar(transitions) > 0))
        return("'G' must have rownames that specify transitions.")

    ## Check that the format of transitions are valid:
    ## For example: "X1 + X2 + ... + Xn -> Y1 + Y2 + ... + Yn"
    ## or
    ## For example: "X1 + X2 + ... + Xn -> propensity -> Y1 + Y2 + ... + Yn"
    ## is expected, where X2, ..., Xn and Y2, ..., Yn are optional.
    transitions <- strsplit(transitions, split = "->", fixed = TRUE)
    if (any(sapply(transitions, length) < 2))
        return("'G' rownames have invalid transitions.")

    ## Check that transitions and S have identical compartments.
    transitions <- unlist(lapply(transitions, function(x) {
        c(x[1], x[length(x)])
    }))
    transitions <- unlist(strsplit(transitions, split = "+", fixed = TRUE))
    transitions <- trimws(transitions)
    transitions <- unique(transitions)
    transitions <- transitions[transitions != "@"]
    transitions <- sub("^[[:digit:]]+[*]", "", transitions)
    if (!all(transitions %in% rownames(object@S)))
        return("'G' and 'S' must have identical compartments.")

    character(0)
}

valid_ldata <- function(object) {
    if (!is.double(object@ldata))
        return("'ldata' matrix must be a double matrix.")
    Nn_ldata <- dim(object@ldata)[2]
    if (Nn_ldata > 0 && !identical(Nn_ldata, dim(object@u0)[2]))
        return("The number of nodes in 'u0' and 'ldata' must match.")

    character(0)
}

valid_gdata <- function(object) {
    if (!is.double(object@gdata))
        return("'gdata' must be a double vector.")

    character(0)
}
