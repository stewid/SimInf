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

#include <string.h>
#include <math.h>

#ifdef SIMINF_OMP
#include <omp.h>
#endif

#include "siminf.h"
#include "events.h"

/**
 * Structure to hold thread specific data/arguments for simulation.
 */
typedef struct siminf_thread_args
{
    int Ni;             /**< Index to first node in thread. */
    int Nn;             /**< Number of nodes in thread. */
    double *data;       /**< Matrix (dsize X Nn). data(:,j) gives a data
                         *   vector for node #j. */
    const int *sd;      /**< Each node can be assigned to a sub-domain. */
    double *sum_t_rate; /**< Vector of length Nn with the sum of propensities
                         *   in every node. */
    double *t_rate;     /**< Transition rate matrix (Nt X Nn) with all
                         *   propensities for state transitions. */
    double *t_time;     /**< Time for next event (transition) in each node. */
    int *state;         /**< Integer vector of length Nn * Nc with state in
                         *   each node. */
    int *individuals;   /**< Vector to store the result of the sampling during
                         *   external events processing. Passed as function
                         *   argument to handle parallellization. */
    int errcode;        /**< The error state of the thread. 0 if ok. */
    const external_events *events; /**< Structure that represents external
                                    * events. */
    int *update_node;   /**< Integer vector of length Nn used to indicate
                         *   nodes for update. */
    gsl_rng *rng;       /**< The random number generator. */
} siminf_thread_args;


/**
 * Initialize transition rate and time to event.
 *
 * Calculate the transition rate for every transition and every
 * node. Store the sum of the transition rates in each node in
 * sum_t_rate. Calculate time to next event (transition) in each
 * node.
 *
 * @param Nn Number of nodes.
 * @param Nc Number of compartments in each node.
 * @param Nt Total number of different transitions.
 * @param dsize Size of data vector sent to propensities.
 * @param state Integer vector of length Nn with state in each node.
 * @param data Double vector (dsize X Nn) with data for each node.
 * @param sd Integer vector of length Nn. Each node can be assigned to
 *        a sub-domain.
 * @param t0 The start time.
 * @param sum_t_rate Double vector of length Nn with the sum of
 *        propensities in every node.
 * @param t_rate Transition rate matrix (Nt X Nn) with all propensities
 *        for state transitions.
 * @param t_time Time for next event (transition) in each node.
 * @param t_fun Vector of function pointers to transition functions.
 * @param rng The random number generator.
 */
static void siminf_init(
    const int Nn, const int Nc, const int Nt, const int dsize,
    const int *state, const double *data, const int *sd,
    const double t0, double *sum_t_rate, double *t_rate,
    double *t_time, const PropensityFun *t_fun, const gsl_rng *rng)
{
    int i, node;

    for (node = 0; node < Nn; node++) {
        sum_t_rate[node] = 0.0;
        for (i = 0; i < Nt; i++) {
            t_rate[node * Nt + i] =
                (*t_fun[i])(
                    &state[node * Nc], t0, &data[node * dsize], sd[node]);

            sum_t_rate[node] += t_rate[node * Nt + i];
        }

        t_time[node] = -log(1.0 - gsl_rng_uniform(rng)) /
            sum_t_rate[node] + t0;
    }
}

/**
 * Epidemiological model
 *
 * Handle internal epidemiological model, continuous-time Markov
 * chain, during one day.
 *
 * G is a sparse matrix dependency graph (Nt X Nt) in
 * compressed column format (CCS). A non-zeros entry in element i of
 * column j indicates that propensity i needs to be recalculated if
 * the transition j occurs.
 *
 * N is a stoichiometry sparse matrix (Nc X Nt) in
 * compressed column format (CCS). Each column corresponds to a
 * transition, and execution of transition j amounts to adding the
 * j'th column to the state vector.
 *
 * @param Nn Number of nodes.
 * @param Nc Number of compartments in each node.
 * @param Nt Total number of different transitions.
 * @param dsize Size of data vector sent to propensities.
 * @param state Integer vector of length Nn with state in each node.
 * @param data Double vector (dsize X Nn) with data for each node.
 * @param sd Integer vector of length Nn. Each node can be assigned to
 *        a sub-domain.
 * @param irG Integer vector where irG[k] is the row of G[k].
 * @param jcG jcG[k], index to data of first non-zero element in row k.
 * @param irN Integer vector where irN[k] is the row of N[k].
 * @param jcN jcN[k], index to data of first non-zero element in row k.
 * @param prN Value of item (i, j) in N.
 * @param sum_t_rate Double vector of length Nn with the sum of
 *        propensities in every node.
 * @param t_rate Transition rate matrix (Nt X Nn) with all propensities
 *        for state transitions.
 * @param t_time Time for next event (transition) in each node.
 * @param next_day Time for next day.
 * @param t_fun Vector of function pointers to transition functions.
 * @param rng The random number generator.
 * @param err The error state of the simulation is saved here.
 *        0 if ok, else error code.
 */
static void siminf_epi_model(
    const int Nn, const int Nc, const int Nt, const int dsize, int *state,
    double *data, const int *sd, const int *irG, const int *jcG,
    const int *irN, const int *jcN, const int *prN, double *sum_t_rate,
    double *t_rate, double *t_time, const double next_day,
    const PropensityFun *t_fun, const gsl_rng *rng, int *err)
{
    int node;

    /* Internal epidemiological model, continuous-time Markov
     * chain. */
    for (node = 0; node < Nn; node++) {
        while (t_time[node] < next_day) {
            double cum, rand, tot_rate, delta = 0.0;
            int i, tr = 0;

            /* a) Determine the transition that did occur (directSSA). */
            cum = t_rate[node * Nt];
            rand = gsl_rng_uniform(rng) * sum_t_rate[node];
            while (tr < Nt && rand > cum) {
                tr++;
                cum += t_rate[node * Nt + tr];
            }

            /* b) Update the state of the node */
            for (i = jcN[tr]; i < jcN[tr + 1]; i++) {
                state[node * Nc + irN[i]] += prN[i];
                if (state[node * Nc + irN[i]] < 0) {
                    *err = SIMINF_ERR_NEGATIVE_STATE;
                    return;
                }
            }

            /* c) Recalculate sum_t_rate[node] using dependency graph. */
            for (i = jcG[tr]; i < jcG[tr + 1]; i++) {
                int j = irG[i];
                double old = t_rate[node * Nt + j];
                delta += (t_rate[node * Nt + j] = (*t_fun[j])(
                              &state[node * Nc], t_time[node],
                              &data[node * dsize], sd[node]))
                    - old;
            }
            sum_t_rate[node] += delta;

            /* d) Compute time to new event for this node. */
            tot_rate = sum_t_rate[node];
            if (tot_rate > 0.0) {
                t_time[node] = -log(1.0 - gsl_rng_uniform(rng)) / tot_rate +
                    t_time[node];
            } else {
                t_time[node] = INFINITY;
            }
        }
    }

    *err = 0;
}

/**
 * Incorporate scheduled external events
 *
 * @param state Integer vector of length Nn * Nc with state in each node.
 * @param event Integer vector of length len with external events.
 * @param time Integer vector of length len with the time for external
          event.
 * @param select Integer vector of length len. Column j in the E
 *        matrix that determines the compartments to sample from.
 * @param node Integer vector of length len. The source node of the
 *        event i.
 * @param dest Integer vector of length len. The dest node of the
 *        event i.
 * @param n Integer vector of length len. The number of individuals
 *        in the external event. n[i] >= 0.
 * @param proportion Double vector of length len. If n[i] equals zero,
 *        then the number of individuals to sample is calculated by
 *        summing the number of individuals in the hidden states
 *        determined by select[i] and multiplying with the proportion.
 *        0 <= p[i] <= 1.
 * @param len Number of scheduled external events.
 * @param index The current index in event list to process. Updated
 *        to the next start index on return.
 * @param irE Array where irE[k] is the row of E[k]
 * @param jcE jcE[k], index to data of first non-zero element in row k
 * @param prE Value of item (i, j) in E.
 * @param individuals The result of the sampling is stored in the
 *        individuals vector. Passed as function argument to handle
 *        parallellization.
 * @param update_node Integer vector of length Nn used to indicate
 *        nodes for update.
 * @param tt The global time.
 * @param rng Random number generator.
 * @param err The error state of processing the events is saved here.
 *        0 if ok, else error code.
 */
static void siminf_process_events(
    const int Nc, const int Nobs, int *state, const int *event,
    const int *time, const int *select, const int *node, const int *dest,
    const int *n, const double *proportion, const int len, int *index,
    const int *irE, const int *jcE, const int *prE, int *individuals,
    int *update_node, const double tt, const gsl_rng *rng, int *err)
{
    int i = *index;
    int e = 0;

    /* Incorporate scheduled external events. */
    while (i < len && tt >= time[i]) {
        e = handle_external_event(
            event[i], irE, jcE, prE, Nc, Nobs, state, node[i], dest[i],
            select[i], n[i], proportion[i], individuals, rng);

        /* Check for error codes. */
        if (e)
            break;

        /* Indicate node and dest node for update */
        update_node[node[i]] = 1;
        update_node[dest[i]] = 1;

        i++;
    }

    *index = i;
    *err = e;
}

/**
 * Update transition rates
 *
 * @param Nn Number of nodes.
 * @param Nc Number of compartments in each node.
 * @param Nt Total number of different transitions.
 * @param dsize Size of data vector sent to propensities.
 * @param state Integer vector of length Nn x Nc with state in each node.
 * @param data Double vector (dsize X Nn) with data for each node.
 * @param sd Integer vector of length Nn. Each node can be assigned to
 *        a sub-domain.
 * @param sum_t_rate Double vector of length Nn with the sum of
 *        propensities in every node.
 * @param t_rate Transition rate matrix (Nt X Nn) with all propensities
 *        for state transitions.
 * @param t_time Time for next event (transition) in each node.
 * @param update_node Integer vector of length Nn to indicate nodes
 *        for update.
 * @param tt The global time.
 * @param t_fun Vector of function pointers to transition functions.
 * @param rng Random number generator.
 */
static void siminf_update(
    const int Nn, const int Nc, const int Nt, const int dsize, int *state,
    double *data, const int *sd, double *sum_t_rate, double *t_rate,
    double *t_time, int *update_node, const double tt,
    const PropensityFun *t_fun, const gsl_rng *rng)
{
    int node;

    for (node = 0; node < Nn; node++) {
        if (update_node[node]) {
            int i = 0;
            double delta = 0.0, old_t_rate = sum_t_rate[node];

            for (; i < Nt; i++) {
                double old = t_rate[node * Nt + i];
                delta += (t_rate[node * Nt + i] =
                          (*t_fun[i])(&state[node * Nc], tt,
                                      &data[node * dsize], sd[node]))
                    - old;
            }
            sum_t_rate[node] += delta;

            if (sum_t_rate[node] > 0.0) {
                if (old_t_rate > 0.0 && !isinf(old_t_rate)) {
                    t_time[node] = old_t_rate / sum_t_rate[node]
                        * (t_time[node] - tt) + tt;
                } else {
                    t_time[node] = -log(1.0 - gsl_rng_uniform(rng))
                        / sum_t_rate[node] + tt;
                }
            } else {
                t_time[node] = INFINITY;
            }

            update_node[node] = 0;
        }
    }
}

/**
 * Single threaded siminf solver
 *
 * G is a sparse matrix dependency graph (Nt X Nt) in
 * compressed column format (CCS). A non-zeros entry in element i of
 * column j indicates that propensity i needs to be recalculated if
 * the event j occurs.
 *
 * N is a stoichiometry sparse matrix (Nc X Nt) in
 * compressed column format (CCS). Each column corresponds to a
 * transition, and execution of transition j amounts to adding the
 * j'th column to the state vector.
 *
 * @param ta Structure (siminf_thread_args) to hold thread specific
 *        data/arguments for simulation.
 * @param irG Integer vector where irG[k] is the row of G[k].
 * @param jcG jcG[k], index to data of first non-zero element in row k.
 * @param irN Integer vector where irN[k] is the row of N[k].
 * @param jcN jcN[k], index to data of first non-zero element in row k.
 * @param prN Value of item (i, j) in N.
 * @param tspan Double vector. Output times. tspan[0] is the start
 *        time and tspan[length(tspan)-1] is the stop time.
 * @param tlen Number of sampling points in time.
 * @param U The output is a matrix U ((Nn * Nc) X length(tspan)).
 *        U(:,j) contains the state of the system at tspan(j).
 * @param Nc Number of compartments in each node.
 * @param Nt Total number of different transitions.
 * @param Nobs Number of observable states.
 * @param dsize Size of data vector sent to propensities.
 * @param irE Integer vector where irE[k] is the row of E[k].
 * @param jcE jcE[k], index to data of first non-zero element in row k.
 * @param prE Value of item (i, j) in E.
 * @param report_level The desired degree of feedback during
 *        simulations. 0, 1, and 2 are currently supported options.
 * @param t_fun Vector of function pointers to transition functions.
 * @param pts_fun Function pointer to callback after each time step
 *        e.g. update infectious pressure.
 * @param progress Function pointer to report progress.
 * @param seed Random number seed.
 * @param update_node Integer vector of length Nn to indicate nodes
 *        for update.
 * @param E1_events external_events structure with E1 events.
 * @param E2_events external_events structure with E2 events.
 * @return 0 if Ok, else error code.
 */
static int siminf_single(
    const int *irG, const int *jcG, const int *irN, const int *jcN,
    const int *prN, const double *tspan, const int tlen, int *U,
    const int *sd, const int Nn, const int Nc, const int Nt,
    const int Nobs, const int dsize, int *state, double *data,
    const int *irE, const int *jcE, const int *prE, int report_level,
    const PropensityFun *t_fun, const PostTimeStepFun pts_fun,
    const ProgressFun progress, unsigned long int seed, int *update_node,
    external_events *E1_events, external_events *E2_events)
{
    double tt = tspan[0];
    long int total_transitions = 0;
    int node, it = 0;
    int errcode = 0;
    double next_day = floor(tspan[0]) + 1.0;
    double *t_rate = NULL, *sum_t_rate = NULL, *t_time = NULL;
    gsl_rng *rng = NULL;
    int *individuals = NULL;

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!rng) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    gsl_rng_set(rng, seed);

    individuals = (int *) malloc(Nc * sizeof(int));
    if (!individuals) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    /* Create transition rate matrix (Nt X Nn) and total rate
     * vector. In t_rate we store all propensities for state
     * transitions, and in sum_t_rate the sum of propensities
     * in every node. */
    t_rate = (double *) malloc(Nt * Nn * sizeof(double));
    if (!t_rate) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    sum_t_rate = (double *) malloc(Nn * sizeof(double));
    if (!sum_t_rate) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    t_time = (double *) malloc(Nn * sizeof(double));
    if (!t_time) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    siminf_init(
        Nn, Nc, Nt, dsize, state, data, sd, tt,
        sum_t_rate, t_rate, t_time, t_fun, rng);

    /* Main loop. */
    for (;;) {
        /* (1) Handle internal epidemiological model, continuous-time
         * Markov chain. */
        siminf_epi_model(
            Nn, Nc, Nt, dsize, state, data, sd, irG,
            jcG, irN, jcN, prN, sum_t_rate, t_rate, t_time,
            next_day, t_fun, rng, &errcode);

        if (errcode)
            break;

        /* (2) Incorporate all scheduled external E1 events. */
        siminf_process_events(
            Nc, Nobs, state, E1_events->event, E1_events->time,
            E1_events->select, E1_events->node, E1_events->dest,
            E1_events->n, E1_events->proportion, E1_events->len,
            &E1_events->index, irE, jcE, prE, individuals,
            update_node, tt, rng, &errcode);

        if (errcode)
            break;

        /* (3) Incorporate all scheduled external E2 events. */
        siminf_process_events(
            Nc, Nobs, state, E2_events->event, E2_events->time,
            E2_events->select, E2_events->node, E2_events->dest,
            E2_events->n, E2_events->proportion, E2_events->len,
            &E2_events->index, irE, jcE, prE, individuals,
            update_node, tt, rng, &errcode);

        if (errcode)
            break;

        /* (4) Incorporate model specific actions after each timestep
         * e.g. update the infectious pressure variable. */
        for (node = 0; node < Nn; node++)
            update_node[node] = pts_fun(&state[node * Nc], node, tt,
                                        &data[node * dsize], sd[node]);

        /* (5) Update transition rates */
        siminf_update(
            Nn, Nc, Nt, dsize, state, data, sd, sum_t_rate,
            t_rate, t_time, update_node, tt, t_fun, rng);

        /* The global time now equals next_day. */
        tt = next_day;

        /* (6) Store solution if tt has passed the next time in *
         * tspan. Report solution up to, but not including tt. */
        if (tt > tspan[it]) {
            for (; it < tlen && tt > tspan[it]; it++) {
                if (report_level)
                    progress(tspan[it], tspan[0], tspan[tlen - 1],
                             total_transitions, report_level);

                memcpy(&U[Nn * Nc * it], state, Nn * Nc * sizeof(int));
            }

            /* If the simulation has reached the final time, exit. */
            if (it >= tlen)
                break;
        }

        next_day += 1.0;
    }


cleanup:
    if (individuals)
        free(individuals);
    if (t_rate)
        free(t_rate);
    if (sum_t_rate)
        free(sum_t_rate);
    if (t_time)
        free(t_time);
    if (rng)
        gsl_rng_free(rng);

    return errcode;
}

/**
 * Initialize and run siminf solver
 *
 * G is a sparse matrix dependency graph (Nt X Nt) in
 * compressed column format (CCS). A non-zeros entry in element i of
 * column j indicates that propensity i needs to be recalculated if
 * the event j occurs.
 *
 * N is a stoichiometry sparse matrix (Nc X Nt) in
 * compressed column format (CCS). Each column corresponds to a
 * transition, and execution of transition j amounts to adding the
 * j'th column to the state vector.
 *
 * @param u0 Initial state vector u0. Integer (Nc X Nn). Gives the
 *        initial number of individuals in each compartment in every
 *        node.
 * @param irG Integer vector where irG[k] is the row of G[k].
 * @param jcG jcG[k], index to data of first non-zero element in row k.
 * @param irN Integer vector where irN[k] is the row of N[k].
 * @param jcN jcN[k], index to data of first non-zero element in row k.
 * @param prN Value of item (i, j) in N.
 * @param tspan Double vector. Output times. tspan[0] is the start
 *        time and tspan[length(tspan)-1] is the stop time.
 * @param tlen Number of sampling points in time.
 * @param U The output is a matrix U ((Nn * Nc) X length(tspan)).
 *        U(:,j) contains the state of the system at tspan(j).
 * @param data Double matrix (dsize X Nn). Generalized data
 *        matrix, data(:,j) gives a data vector for node #j.
 * @param sd Integer vector of length Nn. Each node can be assigned to
 *        a sub-domain.
 * @param Nn Number of nodes.
 * @param Nc Number of compartments in each node.
 * @param Nt Total number of different transitions.
 * @param Nobs Number of observable states.
 * @param dsize Size of data vector sent to propensities.
 * @param irE Integer vector where irE[k] is the row of E[k].
 * @param jcE jcE[k], index to data of first non-zero element in row k.
 * @param prE Value of item (i, j) in E.
 * @param events Structure that represents external events.
 * @param report_level The desired degree of feedback during
 *        simulations. 0, 1, and 2 are currently supported options.
 * @param Nthread Number of threads to use during simulation.
 * @param seed Random number seed.
 * @param t_fun Vector of function pointers to transition functions.
 * @param pts_fun Function pointer to callback after each time step
 *        e.g. update infectious pressure.
 * @param progress Function pointer to report progress.
 * @return 0 if Ok, else error code.
 */
int siminf_run(
    const int *u0, const int *irG, const int *jcG, const int *irN,
    const int *jcN, const int *prN, const double *tspan, const int tlen,
    int *U, double *data, const int *sd, const int Nn, const int Nc,
    const int Nt, const int Nobs, const int dsize, const int *irE,
    const int *jcE, const int *prE, const external_events *events,
    int report_level, int Nthread, unsigned long int seed,
    const PropensityFun *t_fun, const PostTimeStepFun pts_fun,
    const ProgressFun progress)
{
    int err = SIMINF_UNSUPPORTED_PARALLELIZATION;
    int *state = NULL, *update_node = NULL;
    external_events *E1_events = NULL, *E2_events = NULL;

    /* Set state to the initial state. */
    state = (int *) malloc(Nn * Nc * sizeof(int));
    if (!state) {
        err = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    memcpy(state, u0, Nn * Nc * sizeof(int));

    /* Setup vector to keep track of nodes that must be updated due to
     * external events */
    update_node = (int *) calloc(Nn, sizeof(int));
    if (!update_node) {
        err = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    E1_events = calloc(Nthread, sizeof(external_events));
    if (!E1_events) {
        err = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    E2_events = calloc(1, sizeof(external_events));
    if (!E2_events) {
        err = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    err = split_external_events(E1_events, E2_events, events, Nthread);
    if (err)
        goto cleanup;

    if (Nthread == 1) {
        err = siminf_single(
            irG, jcG, irN, jcN, prN, tspan, tlen, U, sd, Nn, Nc, Nt,
            Nobs, dsize, state, data, irE, jcE, prE, report_level,
            t_fun, pts_fun, progress, seed, update_node,
            E1_events, E2_events);
    }

cleanup:
    if (state)
        free(state);
    if (update_node)
        free(update_node);

    if (E1_events) {
        int i;

        for (i = 0; i < Nthread; i++) {
            if (E1_events[i].event)
                free(E1_events[i].event);
            if (E1_events[i].time)
                free(E1_events[i].time);
            if (E1_events[i].select)
                free(E1_events[i].select);
            if (E1_events[i].node)
                free(E1_events[i].node);
            if (E1_events[i].dest)
                free(E1_events[i].dest);
            if (E1_events[i].n)
                free(E1_events[i].n);
            if (E1_events[i].proportion)
                free(E1_events[i].proportion);
        }
        free(E1_events);
    }

    if (E2_events) {
        if (E2_events->event)
            free(E2_events->event);
        if (E2_events->time)
            free(E2_events->time);
        if (E2_events->select)
            free(E2_events->select);
        if (E2_events->node)
            free(E2_events->node);
        if (E2_events->dest)
            free(E2_events->dest);
        if (E2_events->n)
            free(E2_events->n);
        if (E2_events->proportion)
            free(E2_events->proportion);
        free(E2_events);
    }

    return err;
}
