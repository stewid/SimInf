/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
 * Copyright (C) 2015 -- 2019 Stefan Widgren
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

#include <Rdefines.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* This is the maximum number of threads that SimInf will use for
 * OpenMP parallel regions. It is initialised (>= 1) by
 * SimInf_init_threads. */
static int SimInf_max_threads = -1;

/* This is the number of threads that SimInf will use when running a
 * trajectory and depends of both SimInf_max_threads and the model to
 * run. */
static int SimInf_threads = -1;

int SimInf_num_threads()
{
    return SimInf_threads;
}

/* Internal function to specify the number of threads to use in a
 * parallel region. Use all avialable threads if the 'threads'
 * argument is <= 0. */
int SimInf_set_num_threads(int threads)
{
    if (threads <= 0 || threads > SimInf_max_threads)
        threads = SimInf_max_threads;
    SimInf_threads = threads;
    return SimInf_threads;
}

/* Compare and return the minimum value of x, y and the integer value
 * of the environmental variable named 'name' (if it exists). */
#ifdef _OPENMP
static int SimInf_min_env(int x, int y, const char *name)
{
    int z;
    const char *p = getenv(name);

    if (p != NULL && ((z = atoi(p)) < y))
        y = z;
    if (y < x)
        x = y;
    return x;
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
SEXP SimInf_init_threads(SEXP threads)
{
    int old_value = SimInf_max_threads;

#ifdef _OPENMP
    SimInf_max_threads = omp_get_num_procs();

    /* The thread limit can be set with the OMP_THREAD_LIMIT
     * environment variable, for example, CRAN uses
     * OMP_THREAD_LIMIT=2. */
    SimInf_max_threads = SimInf_min_env(
        SimInf_max_threads, omp_get_thread_limit(), "OMP_THREAD_LIMIT");

    /* The maximum number of threads can also be limited with the
     * OMP_NUM_THREADS environment variable. */
    SimInf_max_threads = SimInf_min_env(
        SimInf_max_threads, omp_get_max_threads(), "OMP_NUM_THREADS");

    /* Additionally, the maximum number of threads can be limited with
     * the SIMINF_NUM_THREADS environment variable. */
    SimInf_max_threads = SimInf_min_env(
        SimInf_max_threads, SimInf_max_threads, "SIMINF_NUM_THREADS");

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
    SimInf_max_threads = 1;
#endif

    /* No need to return the number of threads if this function is
     * called from 'R_init_SimInf' during initialisation. */
    return old_value < 1 ? R_NilValue : Rf_ScalarInteger(SimInf_max_threads);
}

/**
 * Is OpenMP available
 */
SEXP SimInf_have_openmp()
{
#ifdef _OPENMP
    return Rf_ScalarLogical(1);
#else
    return Rf_ScalarLogical(0);
#endif
}
