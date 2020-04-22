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

## Raise an error if the error message doesn't match.
check_error <- function(current, target, exact = TRUE) {
    if (isTRUE(exact)) {
        if (!identical(current[[1]]$message, target)) {
            stop("message: ", current[[1]]$message,
                 "\nexpected: ", target)
        }
    } else if (!length(grep(target, current[[1]]$message))) {
        stop("message: ", current[[1]]$message,
             "\nexpected: ", target)
    }

    ## Check that the error message ends with '.'
    if (!length(grep("[.]$", current[[1]]$message))) {
        stop("message: ", current[[1]]$message,
             "\nexpected: ", target)
    }

    invisible(NULL)
}
