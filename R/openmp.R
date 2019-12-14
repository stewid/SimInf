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
##' Set the number of threads to be used in SimInf code that is
##' parallelized with OpenMP (if available). The number of threads is
##' initialized when SimInf is first loaded in the R session using
##' optional envioronment variables (see \sQuote{Details}). It is also
##' possible to specify the number of threads by calling
##' \code{set_num_threads}. If the environment variables that affect
##' the number of threads change, then \code{set_num_threads} must be
##' called again for it to take effect.
##'
##' The \code{omp_get_num_procs()} function is used to determine the
##' number of processors that are available to the device at the time
##' the routine is called. The number of threads is then limited by
##' \code{omp_get_thread_limit()} and the current values of the
##' environmental variables (if set)
##'
##' \itemize{
##'   \item{\code{Sys.getenv("OMP_THREAD_LIMIT")}}
##'   \item{\code{Sys.getenv("OMP_NUM_THREADS")}}
##'   \item{\code{Sys.getenv("SIMINF_NUM_THREADS")}}
##' }
##'
##' Additionally, the maximum number of threads can be controlled by
##' the \code{threads} argument, given that its value is not above any
##' of the limits described above.
##' @param threads integer with maximum number of threads to use in
##'     functions that are parallelized with OpenMP (if
##'     available). Default is NULL, i.e. to use all available
##'     processors and then check for limits in the environment
##'     varibles (see \sQuote{Details}).
##' @return The previous value is returned (invisible).
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
