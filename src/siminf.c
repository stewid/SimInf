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
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <time.h>

#include "siminf.h"
#include "SISe.h"
#include "SISe3.h"

/**
 * Report error
 *
 * @param err The error code.
 */
void siminf_error(int err)
{
    switch (err) {
    case SIMINF_ERR_NEGATIVE_STATE:
        Rf_error("Negative state detected.");
        break;
    case SIMINF_ERR_ALLOC_MEMORY_BUFFER:
        Rf_error("Unable to allocate memory buffer");
        break;
    default:
        Rf_error("Unknown error code.");
    }
}

/**
 * Report progress for siminf solver
 *
 * @param t Current time of simulation.
 * @param t_begin Start time of simulation.
 * @param t_end End time of simulation.
 * @param total_transitions Total number of transition events.
 * @param report_level Level of siminf report during simulation.
 *        Silent if 0, progress if 1 and verbose if 2.
 */
void progress(
    double t,
    const double t_begin,
    const double t_end,
    long int total_transitions,
    int report_level)
{
    switch (report_level) {
    case 1:
        Rprintf("%i%% done.\n",(int)((t - t_begin) / (t_end - t_begin) * 100.0));
        break;
    case 2:
        Rprintf("%i%% done.\n",(int)((t - t_begin) / (t_end - t_begin) * 100.0));
        Rprintf("\t#Transition events = %li\n", total_transitions);
        break;
    default:
        break;
    }
}

/**
 * Get report level
 *
 * @param verbose Level of feedback from simulation
 * @return integer with report level
 */
static int get_report_level(SEXP verbose)
{
    int n;

    if (verbose == R_NilValue)
        Rf_error("verbose must be specified");
    if (LENGTH(verbose) != 1)
        Rf_error("Invalid length of vebose vector");

    if (isInteger(verbose)) {
        if (INTEGER(verbose)[0] == NA_INTEGER)
            Rf_error("Invalid value (NA) for verbose");
        n = INTEGER(verbose)[0];
    } else if (isReal(verbose)) {
        if (REAL(verbose)[0] == NA_REAL)
            Rf_error("Invalid value (NA) for verbose");
        n = (int)(REAL(verbose)[0]);
    } else {
        Rf_error("Invalid type for verbose");
    }

    if (n < 0 || n > 2)
        Rf_error("verbose must be a 0 <= value <= 0");

    return n;
}

/**
 * Look for seed
 *
 * @param seed Random number seed.
 * @return seed
 */
static unsigned long int get_seed(SEXP seed)
{
    if (seed != R_NilValue) {
        switch (LENGTH(seed)) {
        case 0:
            return (unsigned long int)time(NULL);
        case 1:
            if (isInteger(seed)) {
                if (INTEGER(seed)[0] == NA_INTEGER)
                    Rf_error("Invalid value (NA) of seed");
                return (unsigned long int)INTEGER(seed)[0];
            } else if (isReal(seed)) {
                if (REAL(seed)[0] == NA_REAL)
                    Rf_error("Invalid value (NA) of seed");
                return (unsigned long int)REAL(seed)[0];
            }
            Rf_error("Invalid type of seed");
            break;
        default:
            Rf_error("Invalid length of seed");
            break;
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
        Rf_error("Number of threads must be specified");
    if (LENGTH(threads) != 1)
        Rf_error("Invalid length of threads vector");

    if (isInteger(threads)) {
        if (INTEGER(threads)[0] == NA_INTEGER)
            Rf_error("Invalid value (NA) for threads");
        n = INTEGER(threads)[0];
    } else if (isReal(threads)) {
        if (REAL(threads)[0] == NA_REAL)
            Rf_error("Invalid value (NA) for threads");
        n = (int)(REAL(threads)[0]);
    } else {
        Rf_error("Invalid type for threads");
    }

    if (n < 1)
        Rf_error("Number of threads must be a value > 0");

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
get_sparse_matrix_int(size_t **ir, size_t **jc, int **pr, SEXP m)
{
    if (ir) {
        size_t nir = LENGTH(GET_SLOT(m, Rf_install("i")));
        int *xir = INTEGER(GET_SLOT(m, Rf_install("i")));
        int i = 0;
        *ir = (size_t*)malloc(nir * sizeof(size_t));
        if (!(*ir))
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        for (; i < nir; i++)
            (*ir)[i] = (size_t)xir[i];
    }

    if (jc) {
        size_t njc = LENGTH(GET_SLOT(m, Rf_install("p")));
        int *xjc = INTEGER(GET_SLOT(m, Rf_install("p")));
        int i = 0;
        *jc = (size_t*)malloc(njc * sizeof(size_t));
        if (!(*jc))
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        for (; i < njc; i++)
            (*jc)[i] = (size_t)xjc[i];
    }

    if (pr) {
        size_t npr = LENGTH(GET_SLOT(m, Rf_install("x")));
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
 * @param verbose Level of feedback from simulation
 * @param seed Random number seed.
 * @param t_fun Vector of function pointers to transition functions.
 * @param inf_fun Function pointer to update infectious pressure.
 */
int run_internal(
    SEXP result,
    SEXP threads,
    SEXP verbose,
    SEXP seed,
    const PropensityFun *t_fun,
    const InfPressFun inf_fun)
{
    int err = 0, Nobs = 0, report_level, n_threads;
    SEXP events, E, N;
    gsl_rng *rng = NULL;
    size_t *irN = NULL, *jcN = NULL;
    size_t *irG = NULL, *jcG = NULL;
    size_t *irE = NULL, *jcE = NULL;
    int *prE = NULL, *prN = NULL;
    double *data = NULL;
    size_t Nn, Nc, tlen, dsize, Nt;

    /* number of threads */
    n_threads = get_threads(threads);

    /* report level */
    report_level = get_report_level(verbose);

    /* seed */
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, get_seed(seed));

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
    events = GET_SLOT(result, Rf_install("events"));
    E = GET_SLOT(events, Rf_install("E"));
    err = get_sparse_matrix_int(&irE, &jcE, &prE, E);
    if (err)
        goto cleanup;

    /* data */
    dsize = LENGTH(GET_SLOT(result, Rf_install("data")));
    data = (double*)malloc(dsize * sizeof(double));
    if (!data) {
         err = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
         goto cleanup;
    }
    memcpy(data, REAL(GET_SLOT(result, Rf_install("data"))), dsize * sizeof(double));
    dsize = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("data")), R_DimSymbol))[0];

    /* Constants */
    Nn   = INTEGER(GET_SLOT(result, Rf_install("Nn")))[0];
    Nc   = INTEGER(GET_SLOT(N, Rf_install("Dim")))[0];
    Nt   = INTEGER(GET_SLOT(N, Rf_install("Dim")))[1];
    tlen = LENGTH(GET_SLOT(result, Rf_install("tspan")));

    /* Calculate number of observable states from the event
     * matrix. Divide the number of columns by four; the number of
     * event handlers. */
    Nobs = INTEGER(GET_SLOT(E, Rf_install("Dim")))[1] / 4;

    /* Output array (to hold a single trajectory) */
    SET_SLOT(result, Rf_install("U"), allocMatrix(INTSXP, Nn * Nc, tlen));

    /* Core simulation routine. */
    err = siminf_core(
        INTEGER(GET_SLOT(result, Rf_install("u0"))),
        irG, jcG,
        irN, jcN, prN,
        REAL(GET_SLOT(result, Rf_install("tspan"))),
        tlen,
        INTEGER(GET_SLOT(result, Rf_install("U"))),
        data,
        INTEGER(GET_SLOT(result, Rf_install("sd"))),
        Nn, Nc, Nt, Nobs, dsize,
        irE, jcE, prE,
        INTEGER(GET_SLOT(events, Rf_install("ext_event"))),
        INTEGER(GET_SLOT(events, Rf_install("ext_time"))),
        INTEGER(GET_SLOT(events, Rf_install("ext_select"))),
        INTEGER(GET_SLOT(events, Rf_install("ext_node"))),
        INTEGER(GET_SLOT(events, Rf_install("ext_dest"))),
        INTEGER(GET_SLOT(events, Rf_install("ext_n"))),
        REAL(GET_SLOT(events,    Rf_install("ext_p"))),
        INTEGER(GET_SLOT(events, Rf_install("ext_len")))[0],
        report_level, n_threads, rng, t_fun, inf_fun, &progress);

cleanup:
    if (rng)
        gsl_rng_free(rng);
    if (data)
        free(data);
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
    if (prE)
        free(prE);

    return err;
}

/**
 * Run simulation for the SISe model
 *
 * @param model The SISe model
 * @param threads Number of threads
 * @param verbose Level of feedback from simulation
 * @param seed Random number seed.
 * @return S4 class SISe with the simulated trajectory in U
 */
SEXP SISe_run(SEXP model, SEXP threads, SEXP verbose, SEXP seed)
{
    int err = 0;
    SEXP result, class_name;
    PropensityFun t_fun[] = {&SISe_S_to_I, &SISe_I_to_S};

    if (R_NilValue == model || S4SXP != TYPEOF(model))
        Rf_error("Invalid SISe model");

    class_name = getAttrib(model, R_ClassSymbol);
    if (strcmp(CHAR(STRING_ELT(class_name, 0)), "SISe") != 0)
        Rf_error("Invalid SISe model: %s", CHAR(STRING_ELT(class_name, 0)));

    result = PROTECT(duplicate(model));

    err = run_internal(
        result,
        threads,
        verbose,
        seed,
        t_fun,
        &SISe_update_infectious_pressure);

    UNPROTECT(1);

    if (err)
        siminf_error(err);

    return result;
}

/**
 * Run simulation for the SISe3 model
 *
 * @param model The SISe3 model
 * @param threads Number of threads
 * @param verbose Level of feedback from simulation
 * @param seed Random number seed.
 * @return S4 class SISe3 with the simulated trajectory in U
 */
SEXP SISe3_run(SEXP model, SEXP threads, SEXP verbose, SEXP seed)
{
    int err = 0;
    SEXP result, class_name;
    PropensityFun t_fun[] = {&SISe3_S_1_to_I_1,
                             &SISe3_I_1_to_S_1,
                             &SISe3_S_2_to_I_2,
                             &SISe3_I_2_to_S_2,
                             &SISe3_S_3_to_I_3,
                             &SISe3_I_3_to_S_3};

    if (R_NilValue == model || S4SXP != TYPEOF(model))
        Rf_error("Invalid SISe3 model");

    class_name = getAttrib(model, R_ClassSymbol);
    if (strcmp(CHAR(STRING_ELT(class_name, 0)), "SISe3") != 0)
        Rf_error("Invalid SISe3 model: %s", CHAR(STRING_ELT(class_name, 0)));

    result = PROTECT(duplicate(model));

    err = run_internal(
        result,
        threads,
        verbose,
        seed,
        t_fun,
        &SISe3_update_infectious_pressure);

    UNPROTECT(1);

    if (err)
        siminf_error(err);

    return result;
}

static const R_CallMethodDef callMethods[] =
{
    {"SISe_run", (DL_FUNC)&SISe_run, 4},
    {"SISe3_run", (DL_FUNC)&SISe3_run, 4},
    {NULL, NULL, 0}
};

/**
* Register routines to R.
*
* @param info Information about the DLL being loaded
*/
void
R_init_siminf(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
