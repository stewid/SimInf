/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
 * Copyright (C) 2015 -- 2024 Stefan Widgren
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

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "SimInf.h"

/* Declare functions to register */
SEXP SEIR_run(SEXP, SEXP);
SEXP SIR_run(SEXP, SEXP);
SEXP SIS_run(SEXP, SEXP);
SEXP SISe_run(SEXP, SEXP);
SEXP SISe3_run(SEXP, SEXP);
SEXP SISe3_sp_run(SEXP, SEXP);
SEXP SISe_sp_run(SEXP, SEXP);
SEXP SimInf_abc_proposals(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP SimInf_abc_weights(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP SimInf_clean_indiv_events(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP SimInf_distance_matrix(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP SimInf_have_openmp(void);
SEXP SimInf_init_threads(SEXP);
SEXP SimInf_lambertW0(SEXP);
SEXP SimInf_ldata_sp(SEXP, SEXP, SEXP);
SEXP SimInf_split_events(SEXP, SEXP);
SEXP SimInf_systematic_resampling(SEXP);
SEXP SimInf_trajectory(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef callMethods[] =
{
    CALLDEF(SEIR_run, 2),
    CALLDEF(SIR_run, 2),
    CALLDEF(SISe_run, 2),
    CALLDEF(SIS_run, 2),
    CALLDEF(SISe3_run, 2),
    CALLDEF(SISe3_sp_run, 2),
    CALLDEF(SISe_sp_run, 2),
    CALLDEF(SimInf_abc_proposals, 8),
    CALLDEF(SimInf_abc_weights, 7),
    CALLDEF(SimInf_clean_indiv_events, 5),
    CALLDEF(SimInf_distance_matrix, 5),
    CALLDEF(SimInf_have_openmp, 0),
    CALLDEF(SimInf_init_threads, 1),
    CALLDEF(SimInf_lambertW0, 1),
    CALLDEF(SimInf_ldata_sp, 3),
    CALLDEF(SimInf_split_events, 2),
    CALLDEF(SimInf_systematic_resampling, 1),
    CALLDEF(SimInf_trajectory, 10),
    {NULL, NULL, 0}
};

/**
 * Register routines to R.
 *
 * @param info Information about the DLL being loaded
 */
void attribute_visible
R_init_SimInf(
    DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
    R_RegisterCCallable("SimInf", "SimInf_local_spread",
                        (DL_FUNC) &SimInf_local_spread);
    R_RegisterCCallable("SimInf", "SimInf_forward_euler_linear_decay",
                        (DL_FUNC) &SimInf_forward_euler_linear_decay);
    R_RegisterCCallable("SimInf", "SimInf_run",
                        (DL_FUNC) &SimInf_run);
    SimInf_init_threads(R_NilValue);
}
