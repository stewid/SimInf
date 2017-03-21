/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015 - 2017 Stefan Engblom
 *  Copyright (C) 2015 - 2017 Stefan Widgren
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

#ifndef INCLUDE_SIMINF_H
#define INCLUDE_SIMINF_H

#include <R.h>
#include <Rinternals.h>

/* Error constants */
typedef enum {
    SIMINF_ERR_NEGATIVE_STATE       = -1,
    SIMINF_ERR_ALLOC_MEMORY_BUFFER  = -2,
    SIMINF_UNDEFINED_EVENT          = -3,
    SIMINF_INVALID_EDGE_PROBABILITY = -4,
    SIMINF_INVALID_SEED_VALUE       = -5,
    SIMINF_INVALID_THREADS_VALUE    = -6,
    SIMINF_ERR_V_IS_NOT_FINITE      = -7,
    SIMINF_ERR_SAMPLE_SELECT        = -8,
    SIMINF_ERR_INVALID_MODEL        = -9,
    SIMINF_ERR_V_IS_NEGATIVE        = -10,
    SIMINF_ERR_INVALID_RATE         = -11
} SimInf_error_code;

/* Forward declaration of the transition rate function. */
typedef double (*TRFun)(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t);

/* Forward declaration of the post time step callback function. */
typedef int (*PTSFun)(
    double *v_new,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    int node,
    double t);

/* Forward declaration of the function to initiate and run the
 * simulation */
SEXP SimInf_run(
    SEXP model,
    SEXP threads,
    SEXP seed,
    TRFun *tr_fun,
    PTSFun pts_fun);

#endif
