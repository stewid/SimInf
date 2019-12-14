## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2019 Stefan Widgren
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

##' Is OpenMP available
##'
##' @return TRUE if SimInf was built with support for OpenMP, else
##'     FALSE.
##' @noRd
have_openmp <- function() {
    .Call(SimInf_have_openmp)
}


##' Specify the number of threads that SimInf should use
##'
##' Set number of threads to be used in SimInf functions that are
##' parallelized with OpenMP (if available). If OpenMP is available,
##' SimInf uses the 'omp_get_num_procs()' function and the
##' environmental variables 'OMP_THREAD_LIMIT', 'OMP_NUM_THREADS', and
##' 'SIMINF_NUM_THREADS' to determine the maximum number of threads to
##' use in functions that are parallelized. Additionally, the maximum
##' number of threads can be controlled by the 'threads' argument.
##' @param threads integer with maximum number of threads to use in
##'     functions that are parallelized with OpenMP (if
##'     available). Default is NULL, i.e. to use all available
##'     processors and then check for limits in the environment
##'     varibles.
##' @return Integer with the maximum number of threads that will be
##'     used.
##' @export
set_num_threads <- function(threads = NULL) {
    if (!is.null(threads)) {
        if (!is.numeric(threads)) {
            stop("'threads' must be an integer >= 1.", call. = FALSE)
        }

        if (any(length(threads) != 1,
                any(is.na(threads)),
                !all(is_wholenumber(threads)),
                any(threads < 1))) {
            stop("'threads' must be an integer >= 1.", call. = FALSE)
        }

        threads <- as.integer(threads)
    }

    invisible(.Call(SimInf_init_threads, threads))
}
