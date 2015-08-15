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

#ifndef INCLUDE_SIMINF_SOLVER_H
#define INCLUDE_SIMINF_SOLVER_H

#include "siminf.h"

/* Definition of function to initialize and run siminf solver */
int siminf_run(
    const int *u0, const double *v0, const int *irG, const int *jcG,
    const int *irN, const int *jcN, const int *prN, const double *tspan,
    int tlen, int *U, double *V, const double *d0, const int *sd, int Nn,
    int Nc, int Nt, int Nd, int dsize, const int *irE, const int *jcE,
    const int *jcS, const int *prS, int len, const int *event,
    const int *time, const int *node, const int *dest, const int *n,
    const double *proportion, const int *select, const int *shift,
    int Nthread, unsigned long int seed, PropensityFun *t_fun,
    PostTimeStepFun pts_fun);

#endif
