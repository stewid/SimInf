/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
 * Copyright (C) 2015 -- 2026 Stefan Widgren
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

#include "SimInf_internal.h"
#include <R_ext/Visibility.h>
#include <Rdefines.h>

static void
SimInf_raise_error(
    int err)
{
    switch (err) {
    case SIMINF_ERR_NEGATIVE_STATE:
        Rf_error("Negative state detected.");
        break;
    case SIMINF_ERR_ALLOC_MEMORY_BUFFER:       /* #nocov */
        Rf_error("Unable to allocate memory buffer.");  /* #nocov */
        break;
    case SIMINF_UNDEFINED_EVENT:
        Rf_error("Undefined event type.");
        break;
    case SIMINF_ERR_V_IS_NOT_FINITE:
        Rf_error("The continuous state 'v' is not finite.");
        break;
    case SIMINF_ERR_SAMPLE_SELECT:
        Rf_error("Unable to sample individuals for event.");
        break;
    case SIMINF_ERR_INVALID_MODEL:
        Rf_error("Invalid model.");
        break;
    case SIMINF_ERR_V_IS_NEGATIVE:
        Rf_error("The continuous state 'v' is negative.");
        break;
    case SIMINF_ERR_INVALID_RATE:
        Rf_error("Invalid rate detected (non-finite or < 0.0).");
        break;
    case SIMINF_ERR_UNKNOWN_SOLVER:
        Rf_error("Invalid 'solver' value.");
        break;
    case SIMINF_ERR_DEST_OUT_OF_BOUNDS:
        Rf_error("'dest' is out of bounds.");
        break;
    case SIMINF_ERR_NODE_OUT_OF_BOUNDS:
        Rf_error("'node' is out of bounds.");
        break;
    case SIMINF_ERR_EVENTS_N:
        Rf_error("'N' is invalid.");
        break;
    case SIMINF_ERR_EVENT_SHIFT:
        Rf_error("'shift' is invalid.");
        break;
    case SIMINF_ERR_SHIFT_OUT_OF_BOUNDS:
        Rf_error("'shift' is out of bounds.");
        break;
    case SIMINF_ERR_INVALID_PROPORTION:
        Rf_error("Invalid proportion detected (< 0.0 or > 1.0).");
        break;
    case SIMINF_ERR_AEM_REPLICATED_MODEL:
        Rf_error("Cannot run the 'aem' solver on a replicated model.");
        break;
    case SIMINF_ERR_SPARSE_MODEL:
        Rf_error("Invalid sparse matrix in model.");
        break;
    case SIMINF_ERR_INVALID_OPTIONS:
        Rf_error("Invalid options.");
        break;
    case SIMINF_ERR_INVALID_SEED:
        Rf_error("Invalid 'seed' value.");
    case SIMINF_ERR_INVALID_CRN:
        Rf_error("Invalid 'crn' value.");
    case SIMINF_ERR_INVALID_RNG_TYPE:
        Rf_error("Invalid 'rng_type' value.");
    default:                   /* #nocov */
        Rf_error("Unknown error code: %i.", err);       /* #nocov */
        break;
    }
}

/**
 * Get the list element named str (ASCII), or return NULL.
 *
 * From the 'Writing R Extensions' manual.
 */
static SEXP
getListElement(
    SEXP list,
    const char *str)
{
    SEXP elmt = R_NilValue;
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);

    for (int i = 0; i < Rf_length(list); i++) {
        if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            /* ASCII only */
            elmt = VECTOR_ELT(list, i);
            break;
        }
    }

    return elmt;
}

/**
 * Parse solver options from an R list.
 *
 * Extracts solver, crn, rng_type, and seed from the options list and
 * populates the corresponding fields in args. Fields that are not
 * present in the options list retain their zero-initialized default
 * values.  The seed field is parsed last to preserve R's RNG state if
 * any other option validation fails.
 *
 * @param options R list with solver options.
 * @param args Solver arguments structure to populate.
 * @return 0 if OK, else error code.
 */
static int
SimInf_parse_options(
    SEXP options,
    SimInf_solver_args *args)
{
    if (!Rf_isNewList(options))
        return SIMINF_ERR_INVALID_OPTIONS;

    /* Solver */
    SEXP solver = getListElement(options, "solver");
    if (!Rf_isNull(solver)) {
        if (!Rf_isString(solver) ||
            Rf_length(solver) != 1 ||
            STRING_ELT(solver, 0) == NA_STRING)
            return SIMINF_ERR_UNKNOWN_SOLVER;

        const char *s = CHAR(STRING_ELT(solver, 0));
        if (strcmp(s, "ssm") == 0)
            args->solver = 0;
        else if (strcmp(s, "aem") == 0)
            args->solver = 1;
        else
            return SIMINF_ERR_UNKNOWN_SOLVER;
    }

    /* Common random numbers */
    SEXP crn = getListElement(options, "crn");
    if (!Rf_isNull(crn)) {
        if (!Rf_isLogical(crn) || Rf_length(crn) != 1 ||
            LOGICAL(crn)[0] == NA_LOGICAL)
            return SIMINF_ERR_INVALID_CRN;
        args->crn = LOGICAL(crn)[0];
    }

    /* RNG type */
    SEXP rng_type = getListElement(options, "rng_type");
    if (!Rf_isNull(rng_type)) {
        if (!Rf_isString(rng_type) ||
            Rf_length(rng_type) != 1 ||
            STRING_ELT(rng_type, 0) == NA_STRING)
            return SIMINF_ERR_INVALID_RNG_TYPE;

        const char *r = CHAR(STRING_ELT(rng_type, 0));
        if (strcmp(r, "mt19937") == 0)
            args->rng_type = 0;
        else if (strcmp(r, "taus2") == 0)
            args->rng_type = 1;
        else
            return SIMINF_ERR_INVALID_RNG_TYPE;
    }

    /* Seed - parsed last to preserve R's RNG state */
    SEXP seed = getListElement(options, "seed");
    if (Rf_isNull(seed)) {
        GetRNGstate();
        args->seed = (unsigned long int)(unif_rand() * UINT_MAX);
        PutRNGstate();
    } else {
        switch (TYPEOF(seed)) {
        case INTSXP:
            if (Rf_length(seed) != 1 ||
                INTEGER(seed)[0] == NA_INTEGER)
                return SIMINF_ERR_INVALID_SEED;
            args->seed = (unsigned long int)INTEGER(seed)[0];
            break;
        case REALSXP:
            if (Rf_length(seed) != 1 ||
                !R_finite(REAL(seed)[0]))
                return SIMINF_ERR_INVALID_SEED;
            args->seed = (unsigned long int)REAL(seed)[0];
            break;
        default:
            return SIMINF_ERR_INVALID_SEED;
        }
    }

    return 0;
}

/**
 * Initiate and run the simulation
 *
 * @param model The SimInf_model
 * @param options Options for running the numerical solver.
 * @param tr_fun Vector of function pointers to transition rate functions.
 * @param pts_fun Function pointer to callback after each time step
 *        e.g. update infectious pressure.
 */
attribute_hidden SEXP
SimInf_run(
    SEXP model,
    SEXP options,
    TRFun *tr_fun,
    PTSFun pts_fun)
{
    int err = 0, nprotect = 0;
    SEXP result = R_NilValue;
    SEXP ext_events, E, G, N, S, prS;
    SEXP tspan;
    SEXP U, V, U_sparse, V_sparse;
    SimInf_solver_args args = { 0 };

    /* If the model ldata is a 0x0 matrix, i.e. Nld == 0, then use
     * ldata_tmp in the transition rate functions. This is to make
     * &ldata[node * Nld] work in the solvers. The reason for INFINITY
     * is to facilitate for the solvers to detect and raise an error
     * if a model C code uses ldata[0] in the transition rate
     * functions. */
    const double ldata_tmp[1] = { INFINITY };

    /* Parse solver options. */
    err = SimInf_parse_options(options, &args);
    if (err)
        goto cleanup;

    /* Check model validity. */
    if (SimInf_arg_check_model(model)) {
        err = SIMINF_ERR_INVALID_MODEL;
        goto cleanup;
    }

    /* Duplicate model. */
    PROTECT(result = Rf_duplicate(model));
    nprotect++;

    /* Dependency graph */
    PROTECT(G = R_do_slot(result, Rf_install("G")));
    nprotect++;
    args.irG = INTEGER(R_do_slot(G, Rf_install("i")));
    args.jcG = INTEGER(R_do_slot(G, Rf_install("p")));

    /* State change matrix */
    PROTECT(S = R_do_slot(result, Rf_install("S")));
    nprotect++;
    PROTECT(prS = Rf_coerceVector(R_do_slot(S, Rf_install("x")), INTSXP));
    nprotect++;
    args.irS = INTEGER(R_do_slot(S, Rf_install("i")));
    args.jcS = INTEGER(R_do_slot(S, Rf_install("p")));
    args.prS = INTEGER(prS);

    /* tspan */
    PROTECT(tspan = R_do_slot(result, Rf_install("tspan")));
    nprotect++;
    args.tspan = REAL(tspan);

    /* Scheduled events */
    PROTECT(ext_events = R_do_slot(result, Rf_install("events")));
    nprotect++;
    args.len = LENGTH(R_do_slot(ext_events, Rf_install("event")));
    args.event = INTEGER(R_do_slot(ext_events, Rf_install("event")));
    args.time = INTEGER(R_do_slot(ext_events, Rf_install("time")));
    args.node = INTEGER(R_do_slot(ext_events, Rf_install("node")));
    args.dest = INTEGER(R_do_slot(ext_events, Rf_install("dest")));
    args.n = INTEGER(R_do_slot(ext_events, Rf_install("n")));
    args.proportion = REAL(R_do_slot(ext_events, Rf_install("proportion")));
    args.select = INTEGER(R_do_slot(ext_events, Rf_install("select")));
    args.shift = INTEGER(R_do_slot(ext_events, Rf_install("shift")));

    /* Select matrix. */
    PROTECT(E = R_do_slot(ext_events, Rf_install("E")));
    nprotect++;
    args.irE = INTEGER(R_do_slot(E, Rf_install("i")));
    args.jcE = INTEGER(R_do_slot(E, Rf_install("p")));
    args.prE = REAL(R_do_slot(E, Rf_install("x")));

    /* Shift matrix. */
    PROTECT(N = R_do_slot(ext_events, Rf_install("N")));
    nprotect++;
    if (Rf_nrows(N) == INTEGER(R_do_slot(E, Rf_install("Dim")))[0])
        args.N = INTEGER(N);

    /* Constants */
    args.Nrep = INTEGER(R_do_slot(result, Rf_install("replicates")))[0];
    args.Nn =
        INTEGER(R_do_slot
                (R_do_slot(result, Rf_install("u0")),
                 R_DimSymbol))[1] / args.Nrep;
    args.Nc = INTEGER(R_do_slot(S, Rf_install("Dim")))[0];
    args.Nt = INTEGER(R_do_slot(S, Rf_install("Dim")))[1];
    args.Nd =
        INTEGER(R_do_slot(R_do_slot(result, Rf_install("v0")), R_DimSymbol))[0];
    args.Nld =
        INTEGER(R_do_slot
                (R_do_slot(result, Rf_install("ldata")), R_DimSymbol))[0];
    args.tlen = LENGTH(R_do_slot(result, Rf_install("tspan")));

    /* Output array (to hold a single trajectory) */
    PROTECT(U_sparse = R_do_slot(result, Rf_install("U_sparse")));
    nprotect++;
    if (SimInf_sparse
        (U_sparse, (ptrdiff_t) args.Nn * (ptrdiff_t) args.Nc, args.Nrep * args.tlen)) {
        args.irU = INTEGER(R_do_slot(U_sparse, Rf_install("i")));
        args.jcU = INTEGER(R_do_slot(U_sparse, Rf_install("p")));
        args.prU = REAL(R_do_slot(U_sparse, Rf_install("x")));
    } else {
        PROTECT(U =
                Rf_allocMatrix(INTSXP, args.Nn * args.Nc,
                               args.Nrep * args.tlen));
        nprotect++;
        R_do_slot_assign(result, Rf_install("U"), U);
        args.U = INTEGER(R_do_slot(result, Rf_install("U")));
    }

    /* Output array (to hold a single trajectory) */
    PROTECT(V_sparse = R_do_slot(result, Rf_install("V_sparse")));
    nprotect++;
    if (SimInf_sparse
        (V_sparse, (ptrdiff_t) args.Nn * (ptrdiff_t) args.Nd, args.Nrep * args.tlen)) {
        args.irV = INTEGER(R_do_slot(V_sparse, Rf_install("i")));
        args.jcV = INTEGER(R_do_slot(V_sparse, Rf_install("p")));
        args.prV = REAL(R_do_slot(V_sparse, Rf_install("x")));
    } else {
        PROTECT(V =
                Rf_allocMatrix(REALSXP, args.Nn * args.Nd,
                               args.Nrep * args.tlen));
        nprotect++;
        R_do_slot_assign(result, Rf_install("V"), V);
        args.V = REAL(R_do_slot(result, Rf_install("V")));
    }

    /* Initial state. */
    args.u0 = INTEGER(R_do_slot(result, Rf_install("u0")));
    args.v0 = REAL(R_do_slot(result, Rf_install("v0")));

    /* Local data */
    if (args.Nld > 0)
        args.ldata = REAL(R_do_slot(result, Rf_install("ldata")));
    else
        args.ldata = ldata_tmp;

    /* Global data */
    args.gdata = REAL(R_do_slot(result, Rf_install("gdata")));

    /* Function pointers */
    args.tr_fun = tr_fun;
    args.pts_fun = pts_fun;

    /* Specify the number of threads to use. */
    if (args.Nrep > 1) {
        /* Limit threads to the number of replicates in the model. */
        args.Nthread = SimInf_set_num_threads(args.Nrep);
    } else {
        /* Limit threads to the number of nodes in the model. */
        args.Nthread = SimInf_set_num_threads(args.Nn);
    }

    /* Run the simulation solver. */
    if (args.solver == 0) {
        /* ssm */
        if (args.Nrep > 1) {
            if (args.crn == 1)
                err = SimInf_run_solver_mssm_crn(&args);
            else
                err = SimInf_run_solver_mssm(&args);
        } else {
            err = SimInf_run_solver_ssm(&args);
        }
    } else {
        /* aem */
        if (args.Nrep > 1)
            err = SIMINF_ERR_AEM_REPLICATED_MODEL;
        else
            err = SimInf_run_solver_aem(&args);
    }

cleanup:
    if (err)
        SimInf_raise_error(err);

    if (nprotect)
        UNPROTECT(nprotect);

    return result;
}
