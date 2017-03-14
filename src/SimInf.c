/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015 - 2017  Stefan Engblom
 *  Copyright (C) 2015 - 2017  Stefan Widgren
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

#include "siminf_arg.h"
#include "siminf_solver.h"

/**
 * Initiate and run the simulation
 *
 * @param model The siminf_model
 * @param threads Number of threads
 * @param seed Random number seed.
 * @param tr_fun Vector of function pointers to transition rate functions.
 * @param pts_fun Function pointer to callback after each time step
 *        e.g. update infectious pressure.
 */
SEXP SimInf_run(
    SEXP model,
    SEXP threads,
    SEXP seed,
    TRFun *tr_fun,
    PTSFun pts_fun)
{
    int err = 0, n_threads;
    SEXP trajectory, names, result = R_NilValue;
    SEXP ext_events, E, G, N, S, prS;
    SEXP U_sparse, V_sparse;
    int *U = NULL, *irU = NULL, *jcU = NULL;
    double *prU = NULL;
    int *irV = NULL, *jcV = NULL;
    double *V = NULL, *prV = NULL;
    int Nn, Nc, Nt, Nd, Nld, tlen;
    unsigned long int s;

    /* Create a list to hold the result of the simulated trajectory. */
    PROTECT(trajectory = allocVector(VECSXP, 2));
    setAttrib(trajectory, R_NamesSymbol, names = allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("error"));
    SET_STRING_ELT(names, 1, mkChar("model"));

    if (siminf_arg_check_model(model)) {
        err = SIMINF_ERR_INVALID_MODEL;
        goto cleanup;
    }

    /* number of threads */
    err = siminf_get_threads(&n_threads, threads);
    if (err)
        goto cleanup;

    /* seed */
    err =  siminf_get_seed(&s, seed);
    if (err)
        goto cleanup;

    /* Duplicate model and add it to the 'model' item in the
     * trajectory list. */
    SET_VECTOR_ELT(trajectory, 1, result = duplicate(model));

    /* SimInf model */
    G = GET_SLOT(result, Rf_install("G"));
    S = GET_SLOT(result, Rf_install("S"));
    PROTECT(prS = coerceVector(GET_SLOT(S, Rf_install("x")), INTSXP));

    /* Scheduled events */
    ext_events = GET_SLOT(result, Rf_install("events"));
    E = GET_SLOT(ext_events, Rf_install("E"));
    N = GET_SLOT(ext_events, Rf_install("N"));

    /* Constants */
    Nn   = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("u0")), R_DimSymbol))[1];
    Nc   = INTEGER(GET_SLOT(S, Rf_install("Dim")))[0];
    Nt   = INTEGER(GET_SLOT(S, Rf_install("Dim")))[1];
    Nd   = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("v0")), R_DimSymbol))[0];
    Nld  = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("ldata")), R_DimSymbol))[0];
    tlen = LENGTH(GET_SLOT(result, Rf_install("tspan")));

    /* Output array (to hold a single trajectory) */
    U_sparse = GET_SLOT(result, Rf_install("U_sparse"));
    if ((INTEGER(GET_SLOT(U_sparse, Rf_install("Dim")))[0] == (Nn * Nc)) &&
        (INTEGER(GET_SLOT(U_sparse, Rf_install("Dim")))[1] == tlen))
    {
        irU = INTEGER(GET_SLOT(U_sparse, Rf_install("i")));
        jcU = INTEGER(GET_SLOT(U_sparse, Rf_install("p")));
        prU = REAL(GET_SLOT(U_sparse, Rf_install("x")));
    } else {
        SET_SLOT(result, Rf_install("U"), allocMatrix(INTSXP, Nn * Nc, tlen));
        U = INTEGER(GET_SLOT(result, Rf_install("U")));
    }

    /* Output array (to hold a single trajectory) */
    V_sparse = GET_SLOT(result, Rf_install("V_sparse"));
    if ((INTEGER(GET_SLOT(V_sparse, Rf_install("Dim")))[0] == (Nn * Nd)) &&
        (INTEGER(GET_SLOT(V_sparse, Rf_install("Dim")))[1] == tlen))
    {
        irV = INTEGER(GET_SLOT(V_sparse, Rf_install("i")));
        jcV = INTEGER(GET_SLOT(V_sparse, Rf_install("p")));
        prV = REAL(GET_SLOT(V_sparse, Rf_install("x")));
    } else {
        SET_SLOT(result, Rf_install("V"), allocMatrix(REALSXP, Nn * Nd, tlen));
        V = REAL(GET_SLOT(result, Rf_install("V")));
    }

    /* Run simulation solver. */
    err = siminf_run_solver(
        INTEGER(GET_SLOT(result, Rf_install("u0"))),
        REAL(GET_SLOT(result, Rf_install("v0"))),
        INTEGER(GET_SLOT(G, Rf_install("i"))),
        INTEGER(GET_SLOT(G, Rf_install("p"))),
        INTEGER(GET_SLOT(S, Rf_install("i"))),
        INTEGER(GET_SLOT(S, Rf_install("p"))),
        INTEGER(prS),
        REAL(GET_SLOT(result, Rf_install("tspan"))),
        tlen,
        U, irU, jcU, prU,
        V, irV, jcV, prV,
        REAL(GET_SLOT(result, Rf_install("ldata"))),
        REAL(GET_SLOT(result, Rf_install("gdata"))),
        Nn, Nc, Nt, Nd, Nld,
        INTEGER(GET_SLOT(E, Rf_install("i"))),
        INTEGER(GET_SLOT(E, Rf_install("p"))),
        INTEGER(N),
        LENGTH(GET_SLOT(ext_events, Rf_install("event"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("event"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("time"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("node"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("dest"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("n"))),
        REAL(GET_SLOT(ext_events,    Rf_install("proportion"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("select"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("shift"))),
        n_threads, s, tr_fun, pts_fun);

cleanup:
    if (err)
        SET_VECTOR_ELT(trajectory, 0, ScalarInteger(err));

    if (result == R_NilValue)
        UNPROTECT(1);
    else
        UNPROTECT(2);

    return trajectory;
}

/**
 * Is OpenMP available
 */
SEXP siminf_have_openmp()
{
    return Rf_ScalarLogical(
#ifdef _OPENMP
        1
#else
        0
#endif
        );
}
