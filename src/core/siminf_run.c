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
 * Get seed value
 *
 * @param out The seed value.
 * @param seed Random number seed from R.
 * @return 0 if Ok, else error code.
 */
static int get_seed(unsigned long int *out, SEXP seed)
{
    int err = 0;

    if (seed != R_NilValue) {
        if (isInteger(seed) || isReal(seed)) {
            switch (LENGTH(seed)) {
            case 0:
                *out = (unsigned long int)time(NULL);
                break;
            case 1:
                if (isInteger(seed)) {
                    if (INTEGER(seed)[0] == NA_INTEGER)
                        err = SIMINF_INVALID_SEED_VALUE;
                    else
                        *out = (unsigned long int)INTEGER(seed)[0];
                } else if (!R_finite(REAL(seed)[0])) {
                    err = SIMINF_INVALID_SEED_VALUE;
                } else {
                    *out = (unsigned long int)REAL(seed)[0];
                }
                break;
            default:
                err = SIMINF_INVALID_SEED_VALUE;
                break;
            }
        } else {
            err = SIMINF_INVALID_SEED_VALUE;
        }
    } else {
        *out = (unsigned long int)time(NULL);
    }

    return err;
}

/**
 * Get number of threads
 *
 * @param out Number of threads
 * @param threads Number of threads from R
 * @return 0 if Ok, else error code.
 */
static int get_threads(int *out, SEXP threads)
{
    int err = 0;

    if (threads == R_NilValue) {
        *out = 0;
    } else if (isInteger(threads)) {
        if (LENGTH(threads) != 1)
            err = SIMINF_INVALID_THREADS_VALUE;
        else if (INTEGER(threads)[0] == NA_INTEGER)
            err = SIMINF_INVALID_THREADS_VALUE;
        else if (INTEGER(threads)[0] < 0)
            err = SIMINF_INVALID_THREADS_VALUE;
        else
            *out = INTEGER(threads)[0];
    } else if (isReal(threads)) {
        if (LENGTH(threads) != 1)
            err = SIMINF_INVALID_THREADS_VALUE;
        else if (!R_finite(REAL(threads)[0]))
            err = SIMINF_INVALID_THREADS_VALUE;
        else if ((int)(REAL(threads)[0] < 0))
            err = SIMINF_INVALID_THREADS_VALUE;
        else
            *out = (int)(REAL(threads)[0]);
    } else {
        err = SIMINF_INVALID_THREADS_VALUE;
    }

    return err;
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
    SEXP ext_events, E, G, N, S, prN, prS;
    int Nn, Nc, Nt, Nd, Nld, tlen;
    unsigned long int s;

    /* number of threads */
    err = get_threads(&n_threads, threads);
    if (err)
        return err;

    /* seed */
    err =  get_seed(&s, seed);
    if (err)
        return err;

    /* G */
    G = GET_SLOT(result, Rf_install("G"));

    /* N */
    N = GET_SLOT(result, Rf_install("N"));
    PROTECT(prN = coerceVector(GET_SLOT(N, Rf_install("x")), INTSXP));

    /* External events */
    ext_events = GET_SLOT(result, Rf_install("events"));
    E = GET_SLOT(ext_events, Rf_install("E"));
    S = GET_SLOT(ext_events, Rf_install("N"));
    PROTECT(prS = coerceVector(GET_SLOT(S, Rf_install("x")), INTSXP));

    /* Constants */
    Nn   = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("u0")), R_DimSymbol))[1];
    Nc   = INTEGER(GET_SLOT(N, Rf_install("Dim")))[0];
    Nt   = INTEGER(GET_SLOT(N, Rf_install("Dim")))[1];
    Nd   = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("v0")), R_DimSymbol))[0];
    Nld  = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("ldata")), R_DimSymbol))[0];
    tlen = LENGTH(GET_SLOT(result, Rf_install("tspan")));

    /* Output array (to hold a single trajectory) */
    SET_SLOT(result, Rf_install("U"), allocMatrix(INTSXP, Nn * Nc, tlen));
    SET_SLOT(result, Rf_install("V"), allocMatrix(REALSXP, Nn * Nd, tlen));

    /* Run simulation solver. */
    err = siminf_run_solver(
        INTEGER(GET_SLOT(result, Rf_install("u0"))),
        REAL(GET_SLOT(result, Rf_install("v0"))),
        INTEGER(GET_SLOT(G, Rf_install("i"))),
        INTEGER(GET_SLOT(G, Rf_install("p"))),
        INTEGER(GET_SLOT(N, Rf_install("i"))),
        INTEGER(GET_SLOT(N, Rf_install("p"))),
        INTEGER(prN),
        REAL(GET_SLOT(result, Rf_install("tspan"))),
        tlen,
        INTEGER(GET_SLOT(result, Rf_install("U"))),
        REAL(GET_SLOT(result, Rf_install("V"))),
        REAL(GET_SLOT(result, Rf_install("ldata"))),
        REAL(GET_SLOT(result, Rf_install("gdata"))),
        INTEGER(GET_SLOT(result, Rf_install("sd"))),
        Nn, Nc, Nt, Nd, Nld,
        INTEGER(GET_SLOT(E, Rf_install("i"))),
        INTEGER(GET_SLOT(E, Rf_install("p"))),
        INTEGER(GET_SLOT(S, Rf_install("p"))),
        INTEGER(prS),
        LENGTH(GET_SLOT(ext_events, Rf_install("event"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("event"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("time"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("node"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("dest"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("n"))),
        REAL(GET_SLOT(ext_events,    Rf_install("proportion"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("select"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("shift"))),
        n_threads, s, t_fun, pts_fun);

    UNPROTECT(2);

    return err;
}
