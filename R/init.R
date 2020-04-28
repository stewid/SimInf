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

## Utility function to coerce the data.frame to a transposed matrix.
as_t_matrix <- function(x) {
    n_col <- ncol(x)
    n_row <- nrow(x)
    lbl <- colnames(x)
    x <- t(data.matrix(x))
    attributes(x) <- NULL
    dim(x) <- c(n_col, n_row)
    rownames(x) <- lbl
    x
}

init_u0 <- function(u0) {
    if (is.null(u0))
        stop("'u0' is NULL.", call. = FALSE)
    if (is.data.frame(u0))
        u0 <- as_t_matrix(u0)
    if (!all(is.matrix(u0), is.numeric(u0)))
        stop("u0 must be an integer matrix.", call. = FALSE)
    if (!is.integer(u0)) {
        if (!all(is_wholenumber(u0)))
            stop("u0 must be an integer matrix.", call. = FALSE)
        storage.mode(u0) <- "integer"
    }

    u0
}

init_G <- function(G) {
    if (!is.null(G)) {
        if (!is(G, "dgCMatrix"))
            G <- as(G, "dgCMatrix")
    }

    G
}

init_S <- function(S) {
    if (!is.null(S)) {
        if (!is(S, "dgCMatrix"))
            S <- as(S, "dgCMatrix")
    }

    S
}

init_ldata <- function(ldata) {
    if (is.null(ldata))
        ldata <- matrix(numeric(0), nrow = 0, ncol = 0)
    if (is.data.frame(ldata))
        ldata <- as_t_matrix(ldata)
    if (is.integer(ldata))
        storage.mode(ldata) <- "double"

    ldata
}

init_gdata <- function(gdata) {
    if (is.null(gdata))
        gdata <- numeric(0)
    if (is.data.frame(gdata)) {
        if (!identical(nrow(gdata), 1L)) {
            stop("When 'gdata' is a data.frame, it must have one row.",
                 call. = FALSE)
        }
        gdata <- unlist(gdata)
    }

    gdata
}

init_U <- function(U) {
    if (is.null(U)) {
        U <- matrix(integer(0), nrow = 0, ncol = 0)
    } else {
        if (!is.integer(U)) {
            if (!all(is_wholenumber(U)))
                stop("U must be an integer.", call. = FALSE)
            storage.mode(U) <- "integer"
        }

        if (!is.matrix(U)) {
            if (!identical(length(U), 0L))
                stop("U must be equal to 0 x 0 matrix.", call. = FALSE)
            dim(U) <- c(0, 0)
        }
    }

    U
}

init_v0 <- function(v0) {
    if (is.null(v0)) {
        v0 <- matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
        if (is.data.frame(v0))
            v0 <- as_t_matrix(v0)
        if (!all(is.matrix(v0), is.numeric(v0)))
            stop("v0 must be a numeric matrix.", call. = FALSE)

        if (!identical(storage.mode(v0), "double"))
            storage.mode(v0) <- "double"
    }

    v0
}

init_V <- function(V) {
    if (is.null(V)) {
        V <- matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
        if (!is.numeric(V))
            stop("V must be numeric.")

        if (!identical(storage.mode(V), "double"))
            storage.mode(V) <- "double"

        if (!is.matrix(V)) {
            if (!identical(length(V), 0L))
                stop("V must be equal to 0 x 0 matrix.", call. = FALSE)
            dim(V) <- c(0, 0)
        }
    }

    V
}

init_C_code <- function(C_code) {
    if (is.null(C_code))
        C_code <- character(0)

    C_code
}

init_tspan <- function(tspan) {
    if (is(tspan, "Date")) {
        ## Coerce the date vector to a numeric vector as days, where
        ## tspan[1] becomes the day of the year of the first year of
        ## the tspan date vector. The dates are added as names to the
        ## numeric vector.
        t0 <- as.numeric(as.Date(format(tspan[1], "%Y-01-01"))) - 1
        tspan_lbl <- format(tspan, "%Y-%m-%d")
        tspan <- as.numeric(tspan) - t0
        names(tspan) <- tspan_lbl
    } else {
        t0 <- NULL
    }
    storage.mode(tspan) <- "double"

    list(tspan = tspan, t0 = t0)
}
