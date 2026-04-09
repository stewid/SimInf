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
##' Set the number of threads to be used in SimInf code that is
##' parallelized with OpenMP (if available).
##'
##' The number of threads is initialized when SimInf is first loaded
##' in the R session, based on optional environment variables (see
##' \sQuote{Details}). It can also be explicitly set by calling
##' \code{set_num_threads}. If the environment variables affecting the
##' thread count change, \code{set_num_threads} must be called again
##' for the new values to take effect.
##'
##' The function determines the number of available processors using
##' \code{omp_get_num_procs()} and limits the thread count based on
##' \code{omp_get_thread_limit()} and the following environment
##' variables (checked in order of precedence):
##' \itemize{
##'   \item \code{SIMINF_NUM_THREADS}: Specific limit for SimInf.
##'   \item \code{OMP_NUM_THREADS}: General OpenMP limit.
##'   \item \code{OMP_THREAD_LIMIT}: Maximum thread limit.
##' }
##'
##' The \code{threads} argument allows you to override these limits,
##' provided the requested value does not exceed the maximum allowed
##' by the environment or hardware constraints.
##'
##' @param threads Integer specifying the maximum number of threads to
##'     use in OpenMP-parallelized functions. If \code{NULL}
##'     (default), SimInf attempts to use all available processors,
##'     subject to the limits imposed by the environment variables
##'     listed above.
##' @return The previous value of the thread count is returned
##'     invisibly.
##' @export
set_num_threads <- function(threads = NULL) {
    if (!is.null(threads)) {
        if (!is.numeric(threads)) {
            stop("'threads' must be an integer >= 1.", call. = FALSE)
        }

        if (any(length(threads) != 1,
                anyNA(threads),
                !all(is_wholenumber(threads)),
                any(threads < 1))) {
            stop("'threads' must be an integer >= 1.", call. = FALSE)
        }

        threads <- as.integer(threads)
    }

    invisible(.Call(SimInf_init_threads, threads))
}
