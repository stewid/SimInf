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

#ifndef INCLUDE_events_h
#define INCLUDE_events_h

#include <gsl/gsl_rng.h>

/* Definition of function to handle external events. */
typedef int (*ExtEventHandlerFun)(
    const size_t *irE, const size_t *jcE, const int *prE, const size_t Nc,
    const int Nobs, int *state, const int node, const int dest,
    const int select, const int n, const double proportion, int *inividuals,
    const gsl_rng *rng);

int event_exit(
    const size_t *irE, const size_t *jcE, const int *prE, const size_t Nc,
    const int Nobs, int *state, const int node, const int dest,
    const int select, const int n, const double proportion, int *inividuals,
    const gsl_rng *rng);

int event_enter(
    const size_t *irE, const size_t *jcE, const int *prE, const size_t Nc,
    const int Nobs, int *state, const int node, const int dest,
    const int select, const int n, const double proportion, int *inividuals,
    const gsl_rng *rng);

int event_internal_transfer(
    const size_t *irE, const size_t *jcE, const int *prE, const size_t Nc,
    const int Nobs, int *state, const int node, const int dest,
    const int select, const int n, const double proportion, int *inividuals,
    const gsl_rng *rng);

int event_external_transfer(
    const size_t *irE, const size_t *jcE, const int *prE, const size_t Nc,
    const int Nobs, int *state, const int node, const int dest,
    const int select, const int n, const double proportion, int *inividuals,
    const gsl_rng *rng);

#endif
