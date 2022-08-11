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

## Raise an error if the error message doesn't match.
check_error <- function(current, target, exact = TRUE) {
    ## Identify the first error condition.
    i <- min(which(vapply(current, inherits, logical(1), "error")))

    if (isTRUE(exact)) {
        if (!identical(current[[i]]$message, target)) {
            stop("message: ", current[[i]]$message,
                 "\nexpected: ", target)
        }
    } else if (!length(grep(target, current[[i]]$message))) {
        stop("message: ", current[[i]]$message,
             "\nexpected: ", target)
    }

    ## Check that the error message ends with '.'
    if (!length(grep("[.]$", current[[i]]$message))) {
        stop("message: ", current[[i]]$message,
             "\nexpected: ", target)
    }

    invisible(NULL)
}
