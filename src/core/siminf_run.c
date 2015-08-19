/*
 *  siminf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015  Stefan Engblom
 *  Copyright (C) 2015  Stefan Widgren
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <Rdefines.h>

#include <time.h>

#include "siminf_solver.h"

/**
 * Look for seed
 *
 * @param seed Random number seed.
 * @return seed
 */
static unsigned long int get_seed(SEXP seed)
{
    if (seed != R_NilValue) {
        if (isInteger(seed) || isReal(seed)) {
            switch (LENGTH(seed)) {
            case 0:
                return (unsigned long int)time(NULL);
            case 1:
                if (isInteger(seed)) {
                    if (INTEGER(seed)[0] == NA_INTEGER)
                        Rf_error("Invalid value (NA) of seed");
                    return (unsigned long int)INTEGER(seed)[0];
                } else if (isReal(seed)) {
                    if (ISNA(REAL(seed)[0]))
                        Rf_error("Invalid value (NA) of seed");
                    return (unsigned long int)REAL(seed)[0];
                }
                break;
            default:
                Rf_error("Invalid length of seed");
                break;
            }
        } else {
            Rf_error("Invalid type of seed");
        }
    }

    return (unsigned long int)time(NULL);
}

/**
 * Get number of threads
 *
 * @param threads Number of threads
 * @return Integer with number of threads
 */
static int get_threads(SEXP threads)
{
    int n;

    if (threads == R_NilValue)
        return 0;

    if (isInteger(threads)) {
        if (LENGTH(threads) != 1)
            Rf_error("Invalid length of threads vector");
        if (INTEGER(threads)[0] == NA_INTEGER)
            Rf_error("Invalid value (NA) for threads");
        n = INTEGER(threads)[0];
    } else if (isReal(threads)) {
        if (LENGTH(threads) != 1)
            Rf_error("Invalid length of threads vector");
        if (ISNA(REAL(threads)[0]))
            Rf_error("Invalid value (NA) for threads");
        n = (int)(REAL(threads)[0]);
    } else {
        Rf_error("Invalid type for threads");
    }

    if (n < 0)
        Rf_error("Number of threads must be a value >= 0");

    return n;
}

/**
 * Extract data from a sparse integer matrix
 *
 * @param ir If non-null, allocate vector and and copy data from m.
 * @param jc If non-null, allocate vector and and copy data from m.
 * @param pr If non-null, allocate vector and and copy data from m.
 * @param m  The sparse integer matrix.
 * @return 0 on success, else error code.
 */
static int
get_sparse_matrix_int(int **ir, int **jc, int **pr, SEXP m)
{
    if (ir) {
        int nir = LENGTH(GET_SLOT(m, Rf_install("i")));
        int *xir = INTEGER(GET_SLOT(m, Rf_install("i")));
        int i = 0;
        *ir = (int*)malloc(nir * sizeof(int));
        if (!(*ir))
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        for (; i < nir; i++)
            (*ir)[i] = xir[i];
    }

    if (jc) {
        int njc = LENGTH(GET_SLOT(m, Rf_install("p")));
        int *xjc = INTEGER(GET_SLOT(m, Rf_install("p")));
        int i = 0;
        *jc = (int*)malloc(njc * sizeof(int));
        if (!(*jc))
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        for (; i < njc; i++)
            (*jc)[i] = xjc[i];
    }

    if (pr) {
        int npr = LENGTH(GET_SLOT(m, Rf_install("x")));
        double *xpr = REAL(GET_SLOT(m, Rf_install("x")));
        int i = 0;
        *pr = (int*)malloc(npr * sizeof(int));
        if (!(*pr))
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        for (; i < npr; i++)
            (*pr)[i] = (int)xpr[i];
    }

    return 0;
}

/**
 * Initiate and run the simulation
 *
 * @param result The siminf_model
 * @param threads Number of threads
 * @param seed Random number seed.
 * @param t_fun Vector of function pointers to transition functions.
 * @param pts_fun Function pointer to callback after each time step
 *        e.g. update infectious pressure.
 */
int siminf_run(
    SEXP result,
    SEXP threads,
    SEXP seed,
    PropensityFun *t_fun,
    PostTimeStepFun pts_fun)
{
    int err = 0, n_threads;
    SEXP ext_events, E, N, S;
    int *irN = NULL, *jcN = NULL, *prN = NULL;
    int *irG = NULL, *jcG = NULL;
    int *irE = NULL, *jcE = NULL;
    int *jcS = NULL, *prS = NULL;
    int Nn, Nc, Nt, Nd, Nld, elen, tlen;
    unsigned long int s;

    /* number of threads */
    n_threads = get_threads(threads);

    /* seed */
    s = get_seed(seed);

    /* G */
    err = get_sparse_matrix_int(&irG, &jcG, NULL, GET_SLOT(result, Rf_install("G")));
    if (err)
        goto cleanup;

    /* N */
    N = GET_SLOT(result, Rf_install("N"));
    err = get_sparse_matrix_int(&irN, &jcN, &prN, N);
    if (err)
        goto cleanup;

    /* External events */
    ext_events = GET_SLOT(result, Rf_install("events"));
    E = GET_SLOT(ext_events, Rf_install("E"));
    err = get_sparse_matrix_int(&irE, &jcE, NULL, E);
    if (err)
        goto cleanup;
    S = GET_SLOT(ext_events, Rf_install("S"));
    err = get_sparse_matrix_int(NULL, &jcS, &prS, S);
    if (err)
        goto cleanup;

    /* Constants */
    Nn   = INTEGER(GET_SLOT(result, Rf_install("Nn")))[0];
    Nc   = INTEGER(GET_SLOT(N, Rf_install("Dim")))[0];
    Nt   = INTEGER(GET_SLOT(N, Rf_install("Dim")))[1];
    Nd   = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("v0")), R_DimSymbol))[0];
    Nld  = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("ldata")), R_DimSymbol))[0];
    elen = LENGTH(GET_SLOT(ext_events, Rf_install("event"))),
    tlen = LENGTH(GET_SLOT(result, Rf_install("tspan")));

    /* Output array (to hold a single trajectory) */
    SET_SLOT(result, Rf_install("U"), allocMatrix(INTSXP, Nn * Nc, tlen));
    SET_SLOT(result, Rf_install("V"), allocMatrix(REALSXP, Nn * Nd, tlen));

    /* Run simulation solver. */
    err = siminf_run_solver(
        INTEGER(GET_SLOT(result, Rf_install("u0"))),
        REAL(GET_SLOT(result, Rf_install("v0"))),
        irG, jcG,
        irN, jcN, prN,
        REAL(GET_SLOT(result, Rf_install("tspan"))),
        tlen,
        INTEGER(GET_SLOT(result, Rf_install("U"))),
        REAL(GET_SLOT(result, Rf_install("V"))),
        REAL(GET_SLOT(result, Rf_install("ldata"))),
        INTEGER(GET_SLOT(result, Rf_install("sd"))),
        Nn, Nc, Nt, Nd, Nld, irE, jcE, jcS, prS, elen,
        INTEGER(GET_SLOT(ext_events, Rf_install("event"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("time"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("node"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("dest"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("n"))),
        REAL(GET_SLOT(ext_events,    Rf_install("proportion"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("select"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("shift"))),
        n_threads, s, t_fun, pts_fun);

cleanup:
    if (irG)
        free(irG);
    if (jcG)
        free(jcG);
    if (irN)
        free(irN);
    if (jcN)
        free(jcN);
    if (prN)
        free(prN);
    if (irE)
        free(irE);
    if (jcE)
        free(jcE);
    if (jcS)
        free(jcS);
    if (prS)
        free(prS);

    return err;
}
