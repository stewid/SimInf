/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015 - 2016 Stefan Engblom
 *  Copyright (C) 2015 - 2016 Stefan Widgren
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

#include "siminf_error.h"

/* Definition of the propensity function. */
typedef double (*PropensityFun)(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd);

/* Definition of the callback function post one time step. */
typedef int (*PostTimeStepFun)(
    double *v_new,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    int node,
    double t,
    int sd);

/* Definition of function to initiate and run the simulation */
SEXP siminf_run(
    SEXP model,
    SEXP threads,
    SEXP seed,
    PropensityFun *t_fun,
    PostTimeStepFun pts_fun);

#endif
