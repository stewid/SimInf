## siminf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015  Stefan Engblom
## Copyright (C) 2015  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' Check for null arguments
##'
##' Raise an error if any of the arguments are null.
##' @param ... The argumens to check
##' @keywords internal
##' @return invisible(NULL)
check_null_arg <- function(...) {
    arg <- list(...)
    for (i in seq_len(length(arg))) {
        if (is.null(arg[[i]]))
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "' is missing"))
    }

    invisible(NULL)
}

##' Check numeric arguments
##'
##' Raise an error if any of the arguments are non-numeric.
##' @param ... The argumens to check
##' @keywords internal
##' @return invisible(NULL)
check_numeric_arg <- function(...) {
    arg <- list(...)
    for (i in seq_len(length(arg))) {
        if (!is.numeric(arg[[i]]))
            stop(paste0("'",
                        match.call(expand.dots = FALSE)$'...'[i],
                        "' must be numeric"))
    }

    invisible(NULL)
}
