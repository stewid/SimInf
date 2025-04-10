## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2025 Stefan Widgren
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

init_error <- function(x, storage_mode) {
    article <- ifelse(identical(storage_mode, "integer"), "an", "a")

    stop(paste0("'", x, "' must be ", article, " ", storage_mode, " matrix."),
         call. = FALSE)
}

init_x0 <- function(x, storage_mode = c("integer", "double"), null_ok = FALSE) {
    storage_mode <- match.arg(storage_mode)

    if (is.null(x)) {
        if (!isTRUE(null_ok))
            init_error(match.call()$x, storage_mode)
        x <- matrix(numeric(0), nrow = 0, ncol = 0)
    }

    if (is.data.frame(x))
        x <- as_t_matrix(x)

    if (!all(is.matrix(x), is.numeric(x)))
        init_error(match.call()$x, storage_mode)

    if (identical(storage_mode, "integer") &&
        !is.integer(x) &&
        !all(is_wholenumber(x))) {
        init_error(match.call()$x, storage_mode)
    }

    if (!identical(storage.mode(x), storage_mode))
        storage.mode(x) <- storage_mode

    x
}

init_sparse_matrix <- function(x) {
    if (!is.null(x) && !methods::is(x, "dgCMatrix")) {
        x <- Matrix::Matrix(x)
        x <- methods::as(x, "dMatrix")
        x <- methods::as(x, "generalMatrix")
        x <- methods::as(x, "CsparseMatrix")
    }

    x
}

init_data_matrix <- function(x) {
    if (is.null(x))
        x <- matrix(numeric(0), nrow = 0, ncol = 0)
    if (is.vector(x, mode = "numeric"))
        x <- as.matrix(x)
    if (is.data.frame(x))
        x <- as_t_matrix(x)
    if (is.integer(x))
        storage.mode(x) <- "double"

    x
}

init_data_vector <- function(x) {
    if (is.null(x))
        x <- numeric(0)
    if (is.data.frame(x)) {
        if (!identical(nrow(x), 1L)) {
            stop(paste0("When '",
                        match.call()$x,
                        "' is a data.frame, it must have one row."),
                 call. = FALSE)
        }
        x <- unlist(x)
    }

    if (is.integer(x))
        storage.mode(x) <- "double"

    x
}

init_output_matrix <- function(x, storage_mode = c("integer", "double")) {
    storage_mode <- match.arg(storage_mode)

    if (is.null(x))
        x <- matrix(numeric(0), nrow = 0, ncol = 0)

    if (!is.numeric(x))
        init_error(match.call()$x, storage_mode)

    if (identical(storage_mode, "integer") &&
        !is.integer(x) &&
        !all(is_wholenumber(x))) {
        init_error(match.call()$x, storage_mode)
    }

    if (!identical(storage.mode(x), storage_mode))
        storage.mode(x) <- storage_mode

    if (!is.matrix(x)) {
        if (!identical(length(x), 0L)) {
            stop(paste0("'", match.call()$x,
                        "' must be equal to a 0 x 0 matrix."),
                 call. = FALSE)
        }
        dim(x) <- c(0, 0)
    }

    x
}

init_C_code <- function(C_code) {
    if (is.null(C_code))
        C_code <- character(0)

    C_code
}

init_tspan <- function(tspan) {
    if (methods::is(tspan, "Date")) {
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
