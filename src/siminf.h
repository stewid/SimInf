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

/* Error constants */
#define SIMINF_ERR_NEGATIVE_STATE          1
#define SIMINF_ERR_ALLOC_MEMORY_BUFFER     2
#define SIMINF_UNSUPPORTED_PARALLELIZATION 3

/* Definition of the propensity function. */
typedef double (*PropensityFun)(const int *x, double t, const double *data, int sd);

/* Definition of the propensity function for infectious pressure. */
typedef int (*InfPressFun)(const int *x, int src, double t, double *data);

/* Definition of the callback function for reporting progress. */
typedef void (*ProgressFun)(double t, const double t_begin, const double t_end,
                            long int total_transitions, int report_level);

int siminf_core(
    const int *u0, const size_t *irG, const size_t *jcG, const size_t *irN,
    const size_t *jcN, const int *prN, const double *tspan, const size_t tlen,
    int *U, double *data, const int *sd, const size_t Nn,
    const size_t Nc, const size_t Nt, const int Nobs, const size_t dsize,
    const size_t *irE, const size_t *jcE, const int *prE, const int *ext_event,
    const int *ext_time, const int *ext_select, const int *ext_node,
    const int *ext_dest, const int *ext_n, const double *ext_p, int ext_len,
    int report_level, int Nthreads, const gsl_rng *rng,
    const PropensityFun *t_fun, const InfPressFun inf_fun,
    const ProgressFun progress);

#endif
