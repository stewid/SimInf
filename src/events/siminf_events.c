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

#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "siminf_arg.h"
#include "siminf_error.h"
#include "siminf_events.h"
#include "siminf_vec.h"

/**
 * Create a directed random graph G_{n,p} without loops
 *
 * Reference: Batagelj, V. & Brandes, U. Efficient generation of large
 * random networks Physical Review E, APS, 2005, 71, 036113
 * @param n Number of nodes
 * @param p Probability of edge
 * @param from The vector to add the 'from' nodes to
 * @param to The vector to add the 'to' nodes to
 * @param rng The random number generator
 * @return 0 if Ok, else error code.
*/
static int siminf_gnp(
    int n,
    double p,
    struct siminf_vec *from,
    struct siminf_vec *to,
    gsl_rng *rng)
{
    int w = -1;
    int v = 0;
    double lp;

    if (p <= 0.0 || p >= 1.0)
        return SIMINF_INVALID_EDGE_PROBABILITY;

    lp = log(1.0 - p);

    while (v < n) {
        double lr = log(1.0 - gsl_rng_uniform(rng));
        w = w + 1 + lr / lp;

        /* Protect against loops */
        if (v == w)
            w += 1;

        while  (w >= n && v < n) {
            w = w - v;
            v = v + 1;

            /* Protect against loops */
            if (v == w)
                w += 1;
        }

        if (v < n) {
            int err;

            err = siminf_vec_push_back(from, v);
            if (err)
                return err;

            err = siminf_vec_push_back(to, w);
            if (err)
                return err;
        }
    }

    return 0;
}

/* Dynamic vectors for events */
struct siminf_events {
    struct siminf_vec event;
    struct siminf_vec time;
    struct siminf_vec node;
    struct siminf_vec dest;
    struct siminf_vec n;
    struct siminf_vec select;
    struct siminf_vec shift;
} siminf_events;

#define SIMINF_EVENTS_INIT {SIMINF_VEC_INIT, \
                            SIMINF_VEC_INIT, \
                            SIMINF_VEC_INIT, \
                            SIMINF_VEC_INIT, \
                            SIMINF_VEC_INIT, \
                            SIMINF_VEC_INIT, \
                            SIMINF_VEC_INIT}

/**
 * Free memory for events
 *
 * @param events The events to free memory for
 * @return void
*/
static void siminf_events_free(struct siminf_events *events)
{
    siminf_vec_free(&(events->event));
    siminf_vec_free(&(events->time));
    siminf_vec_free(&(events->node));
    siminf_vec_free(&(events->dest));
    siminf_vec_free(&(events->n));
    siminf_vec_free(&(events->select));
    siminf_vec_free(&(events->shift));
}

/**
 * Reserves memory for events
 *
 * @param events The events to reserve memory for
 * @param capacity The new capacity of the events
 * @return 0 if Ok, else error code.
*/
static int siminf_events_reserve(struct siminf_events *events, size_t capacity)
{
    int err;

    err = siminf_vec_reserve(&(events->event), capacity);
    if (err)
        return err;
    err = siminf_vec_reserve(&(events->time), capacity);
    if (err)
        return err;
    err = siminf_vec_reserve(&(events->node), capacity);
    if (err)
        return err;
    err = siminf_vec_reserve(&(events->dest), capacity);
    if (err)
        return err;
    err = siminf_vec_reserve(&(events->n), capacity);
    if (err)
        return err;
    err = siminf_vec_reserve(&(events->select), capacity);
    if (err)
        return err;
    err = siminf_vec_reserve(&(events->shift), capacity);
    if (err)
        return err;

    return 0;
}

/**
 * Create random external transfer events
 *
 * @param nodes Number of nodes
 * @param p_edge Probability of edge
 * @param mu The mean number individuals in a transfer event
 * @param t The time of the transfer event
 * @param events The new transfer events are added to events
 * @param rng The random number generator
 * @return 0 if Ok, else error code.
*/
static int siminf_external_transfer_events(
    int nodes,
    double p_edge,
    double mu,
    int t,
    struct siminf_events *events,
    gsl_rng *rng)
{
    int err = siminf_gnp(nodes, p_edge, &(events->node), &(events->dest), rng);
    if (err)
        return err;

    while (events->time.size < events->node.size) {
        int n;

        err = siminf_vec_push_back(&(events->event), 3);
        if (err)
            return err;

        err = siminf_vec_push_back(&(events->time), t);
        if (err)
            return err;

        n = gsl_ran_poisson(rng, mu);
        if (!n)
            n = 1;
        err = siminf_vec_push_back(&(events->n), n);
        if (err)
            return err;

        err = siminf_vec_push_back(&(events->select), 0);
        if (err)
            return err;

        err = siminf_vec_push_back(&(events->shift), -1);
        if (err)
            return err;
    }

    return 0;
}

/**
 * Create random external events
 *
 * @param nodes Number of nodes
 * @param days Number of days
 * @param p_edge Vector of length 'days' with probabilities of
 * edges. One value for each day.
 * @param mu Vector of length 'days' with the mean number individuals
 * in a transfer event. One value for each day.
 * @param seed Random number seed.
 * @return A named list of vectors
*/
SEXP siminf_external_events(
    SEXP nodes,
    SEXP days,
    SEXP p_edge,
    SEXP mu,
    SEXP seed)
{
    int err;
    SEXP result = R_NilValue;
    SEXP names = R_NilValue;
    SEXP item;
    size_t i;
    struct siminf_events events = SIMINF_EVENTS_INIT;
    gsl_rng *rng = NULL;
    size_t capacity = 10;

    /* Check arguments */
    if (siminf_arg_check_integer(nodes))
        Rf_error("Invalid 'nodes' argument");
    if (siminf_arg_check_integer(days))
        Rf_error("Invalid 'days' argument");
    if (siminf_arg_check_real_vec(p_edge, INTEGER(days)[0]))
        Rf_error("Invalid 'p_edge' argument");
    if (siminf_arg_check_real_vec(mu, INTEGER(days)[0]))
        Rf_error("Invalid 'mu' argument");
    if ((seed != R_NilValue) && siminf_arg_check_integer(seed))
        Rf_error("Invalid 'seed' argument");

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (seed == R_NilValue)
        gsl_rng_set(rng, (unsigned long int)time(NULL));
    else
        gsl_rng_set(rng, (unsigned long int)INTEGER(seed)[0]);

    err = siminf_events_reserve(&events, capacity);
    if (err)
        goto cleanup;

    for (i = 0; i < INTEGER(days)[0]; i++) {
        err = siminf_external_transfer_events(
            INTEGER(nodes)[0],
            REAL(p_edge)[i],
            REAL(mu)[i],
            i,
            &events,
            rng);
        if (err)
            goto cleanup;
    }

    /* Copy the events to a named list of vectors. */
    PROTECT(result = allocVector(VECSXP, 7));
    setAttrib(result, R_NamesSymbol, names = allocVector(STRSXP, 7));

    SET_VECTOR_ELT(result, 0, item = allocVector(INTSXP, events.event.size));
    memcpy(INTEGER(item), events.event.buf, events.event.size * sizeof(int));
    SET_STRING_ELT(names, 0, mkChar("event"));

    SET_VECTOR_ELT(result, 1, item = allocVector(INTSXP, events.time.size));
    memcpy(INTEGER(item), events.time.buf, events.time.size * sizeof(int));
    SET_STRING_ELT(names, 1, mkChar("time"));

    SET_VECTOR_ELT(result, 2, item = allocVector(INTSXP, events.node.size));
    memcpy(INTEGER(item), events.node.buf, events.node.size * sizeof(int));
    SET_STRING_ELT(names, 2, mkChar("node"));

    SET_VECTOR_ELT(result, 3, item = allocVector(INTSXP, events.dest.size));
    memcpy(INTEGER(item), events.dest.buf, events.dest.size * sizeof(int));
    SET_STRING_ELT(names, 3, mkChar("dest"));

    SET_VECTOR_ELT(result, 4, item = allocVector(INTSXP, events.n.size));
    memcpy(INTEGER(item), events.n.buf, events.n.size * sizeof(int));
    SET_STRING_ELT(names, 4, mkChar("n"));

    SET_VECTOR_ELT(result, 5, item = allocVector(INTSXP, events.select.size));
    memcpy(INTEGER(item), events.select.buf, events.select.size * sizeof(int));
    SET_STRING_ELT(names, 5, mkChar("select"));

    SET_VECTOR_ELT(result, 6, item = allocVector(INTSXP, events.shift.size));
    memcpy(INTEGER(item), events.shift.buf, events.shift.size * sizeof(int));
    SET_STRING_ELT(names, 6, mkChar("shift"));

cleanup:
    if (rng)
        gsl_rng_free(rng);

    siminf_events_free(&events);

    if (result != R_NilValue)
        UNPROTECT(1);

    if (err)
        siminf_error(err);

    return result;
}
