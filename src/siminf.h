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

#ifndef INCLUDE_siminf_h
#define INCLUDE_siminf_h

#include <stddef.h>
#include <gsl/gsl_rng.h>

#include "events.h"

/* Error constants */
#define SIMINF_ERR_NEGATIVE_STATE          1
#define SIMINF_ERR_ALLOC_MEMORY_BUFFER     2
#define SIMINF_UNSUPPORTED_PARALLELIZATION 3
#define SIMINF_UNDEFINED_EVENT             4

/* Definition of the propensity function. */
typedef double (*PropensityFun)(const int *x, double t, const double *data,
                                int sd);

/* Definition of the callback function post one time step. */
typedef int (*PostTimeStepFun)(const int *x, int src, double t, double *data,
                               int sd);

/* Definition of the callback function for reporting progress. */
typedef void (*ProgressFun)(double t, const double t_begin, const double t_end,
                            long int total_transitions, int report_level);

int siminf_run(
    const int *u0, const int *irG, const int *jcG, const int *irN,
    const int *jcN, const int *prN, const double *tspan, const int tlen,
    int *U, double *data, const int *sd, const int Nn,
    const int Nc, const int Nt, const int Nobs, const int dsize,
    const int *irE, const int *jcE, const int *prE,
    const external_events *events,
    int report_level, int Nthreads, unsigned long int seed,
    const PropensityFun *t_fun, const PostTimeStepFun pts_fun,
    const ProgressFun progress);

#endif
