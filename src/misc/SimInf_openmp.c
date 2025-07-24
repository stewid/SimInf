/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
 * Copyright (C) 2015 -- 2025 Stefan Widgren
 *
 * SimInf is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SimInf is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "SimInf.h"
#ifdef _OPENMP
#  include <omp.h>
#endif
#include <R_ext/Visibility.h>
#include <Rinternals.h>
#include <stdlib.h>

/* This is the maximum number of threads that SimInf will use for
 * OpenMP parallel regions. It is initialised (>= 1) by
 * SimInf_init_threads. */
static int SimInf_max_threads = -1;

/* This is the number of threads that SimInf will use when running a
 * trajectory and depends of both SimInf_max_threads and the model to
 * run. */
static int SimInf_threads = -1;

attribute_hidden
int
SimInf_num_threads(
    void)
{
    return SimInf_threads;
}

/* Internal function to specify the number of threads to use in a
 * parallel region. Use all avialable threads if the 'threads'
 * argument is <= 0. */
attribute_hidden
int
SimInf_set_num_threads(
    int threads)
{
    if (threads <= 0 || threads > SimInf_max_threads)
        threads = SimInf_max_threads;
    SimInf_threads = threads;
    return SimInf_threads;
}

/* Get the value of the environmental variable 'SIMINF_NUM_THREADS'
 * (if it exists and is greater than 0). */
#ifdef _OPENMP
static int SimInf_get_max_threads(void)
{
    const char *p = getenv("SIMINF_NUM_THREADS");

    if (p != NULL) {
        int value = atoi(p);
        if (value > 0)
            return value;
    }

    return INT_MAX;
}
#endif

/* Initialise the maximum number of threads to use in a parallel
 * regions. This function is called from 'R_init_SimInf' but can also
 * be called from R to determine the maximum number of threads to use
 * in a parallel region. It uses information from the
 * omp_get_num_procs() function and the environmental variables
 * 'OMP_THREAD_LIMIT', 'OMP_NUM_THREADS', and 'SIMINF_NUM_THREADS' to
 * find the number of threads. Additionally, it can be controlled by
 * the 'threads' argument when called from 'R'. If called from R, it
 * returns the old value of the maximum number of threads used. */
attribute_hidden
SEXP
SimInf_init_threads(
    SEXP threads)
{
    int old_value = SimInf_max_threads;

#ifdef _OPENMP
    int thread_limit;
    SimInf_max_threads = omp_get_num_procs();

    /* The thread limit can be set with the OMP_THREAD_LIMIT
     * environment variable, for example, CRAN uses
     * OMP_THREAD_LIMIT=2. */
    if ((thread_limit = omp_get_thread_limit()) < SimInf_max_threads)
        SimInf_max_threads = thread_limit;

    /* The maximum number of threads can also be limited with the
     * OMP_NUM_THREADS environment variable. */
    if ((thread_limit = omp_get_max_threads()) < SimInf_max_threads)
        SimInf_max_threads = thread_limit;

    /* Additionally, the maximum number of threads can be limited with
     * the SIMINF_NUM_THREADS environment variable. */
    if ((thread_limit = SimInf_get_max_threads()) < SimInf_max_threads)
        SimInf_max_threads = thread_limit;

    if (Rf_isInteger(threads) &&
        LENGTH(threads) == 1 &&
        INTEGER(threads)[0] != NA_INTEGER &&
        INTEGER(threads)[0] < SimInf_max_threads)
    {
        SimInf_max_threads = INTEGER(threads)[0];
    }

    /* Make sure to have at least one available thread. */
    if (SimInf_max_threads < 1)
        SimInf_max_threads = 1;
#else
    SIMINF_UNUSED(threads);
    SimInf_max_threads = 1;
#endif

    /* No need to return the number of threads if this function is
     * called from 'R_init_SimInf' during initialisation. */
    return old_value < 1 ? R_NilValue : Rf_ScalarInteger(old_value);
}

/**
 * Is OpenMP available
 */
attribute_hidden
SEXP
SimInf_have_openmp(
    void)
{
#ifdef _OPENMP
    return Rf_ScalarLogical(1);
#else
    return Rf_ScalarLogical(0);
#endif
}
