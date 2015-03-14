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
 * @param state Integer vector of length Nn with state in each node.
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
 * Post timestep
 *
 * Incorporate model specific actions after each timestep e.g. update
 * the infectious pressure variable.
 *
 * @param Nn Number of nodes.
 * @param Nc Number of compartments in each node.
 * @param Nt Total number of different transitions.
 * @param dsize Size of data vector sent to propensities.
 * @param state Integer vector of length Nn with state in each node.
 * @param data Double vector (dsize X Nn) with data for each node.
 * @param sd Integer vector of length Nn. Each node can be assigned to
 *        a sub-domain.
 * @param sum_t_rate Double vector of length Nn with the sum of
 *        propensities in every node.
 * @param t_rate Transition rate matrix (Nt X Nn) with all propensities
 *        for state transitions.
 * @param t_time Time for next event (transition) in each node.
 * @param update_node Integer vector of length Nn used to indicate
 *        nodes for update.
 * @param tt The global time.
 * @param t_fun Vector of function pointers to transition functions.
 * @param pts_fun Function pointer to callback after each time step
 *        e.g. update infectious pressure.
 * @param rng Random number generator.
 */
static void siminf_post_timestep(
    const int Nn, const int Nc, const int Nt, const int dsize, int *state,
    double *data, const int *sd, double *sum_t_rate, double *t_rate,
    double *t_time, int *update_node, const double tt,
    const PropensityFun *t_fun, const PostTimeStepFun pts_fun,
    const gsl_rng *rng)
{
    int node;

    for (node = 0; node < Nn; node++) {
        if (pts_fun(
                &state[node * Nc], node, tt, &data[node * dsize], sd[node]) ||
            update_node[node])
        {
            int i = 0;
            double delta = 0.0, old_t_rate = sum_t_rate[node];

            /* compute new transition rate only for transitions
             * dependent on infectious pressure */
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
 * Core siminf solver
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
 * @param Nthread Number of threads to use during simulation. Always 1
 *        for this solver.
 * @param rng The random number generator.
 * @param t_fun Vector of function pointers to transition functions.
 * @param pts_fun Function pointer to callback after each time step
 *        e.g. update infectious pressure.
 * @param progress Function pointer to report progress.
 */

static int siminf_core_single(
    const int *u0, const int *irG, const int *jcG, const int *irN,
    const int *jcN, const int *prN, const double *tspan, const int tlen,
    int *U, double *data, const int *sd, const int Nn,
    const int Nc, const int Nt, const int Nobs, const int dsize,
    const int *irE, const int *jcE, const int *prE,
    const external_events *events,
    int report_level, int Nthreads, const gsl_rng *rng,
    const PropensityFun *t_fun, const PostTimeStepFun pts_fun,
    const ProgressFun progress)
{
    double tt = tspan[0];
    double *sum_t_rate = NULL, *t_rate = NULL, *t_time = NULL;
    int *xx = NULL;
    int *individuals = NULL;
    long int total_transitions = 0;
    int node, errcode = 0;
    int it = 0;
    const int Ndofs = Nn * Nc;

    /* Variables to handle external events */
    int ext_i                    = 0;
    double next_day = floor(tspan[0]) + 1.0;
    int *update_node = NULL;

    individuals = malloc(Nc * sizeof(int));
    if (!individuals) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    /* Setup vector to keep track of nodes that must be updated due to
     * external events */
    update_node = (int *) calloc(Nn, sizeof(int));
    if (!update_node) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    /* Set xx to the initial state. */
    xx = (int *) malloc(Ndofs * sizeof(int));
    if (!xx) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    memcpy(xx, u0, Ndofs * sizeof(int));

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

    /* Calculate the transition rate for every transition and every
     * node. Store the sum of the transition rates in each node in
     * sum_t_rate. Calculate time to next event (transition) in each
     * node. */
    for (node = 0; node < Nn; node++) {
        int i;

        sum_t_rate[node] = 0.0;
        for (i = 0; i < Nt; i++) {
            t_rate[node * Nt + i] = (*t_fun[i])(&xx[node * Nc],
                                                tt,
                                                &data[node * dsize],
                                                sd[node]);

            sum_t_rate[node] += t_rate[node * Nt + i];
        }

        t_time[node] = -log(1.0 - gsl_rng_uniform(rng)) / sum_t_rate[node] +
            tspan[0];
    }

    /* Main loop. */
    for (;;) {
        /* (1) Handle internal epidemiological model, continuous-time
         * Markov chain. */
        siminf_epi_model(
            Nn, Nc, Nt, dsize, xx, data, sd, irG, jcG, irN, jcN, prN,
            sum_t_rate, t_rate, t_time, next_day, t_fun, rng, &errcode);

        if (errcode)
            break;

        /* (2) Incorporate all scheduled external events. */
        siminf_process_events(
            Nc, Nobs, xx, events->event, events->time, events->select,
            events->node, events->dest, events->n, events->proportion,
            events->len, &ext_i, irE, jcE, prE, individuals,
            update_node, tt, rng, &errcode);

        if (errcode)
            break;

        /* (3) Incorporate model specific actions after each timestep
         * e.g. update the infectious pressure variable. */
        siminf_post_timestep(
            Nn, Nc, Nt, dsize, xx, data, sd, sum_t_rate, t_rate,
            t_time, update_node, tt, t_fun, pts_fun, rng);

        /* The global time now equals next_day. */
        tt = next_day;

        /* Store solution if tt has passed the next time in *
         * tspan. Report solution up to, but not including tt. */
        if (tt > tspan[it]) {
            for (; it < tlen && tt > tspan[it]; it++) {
                if (report_level)
                    progress(tspan[it], tspan[0], tspan[tlen - 1],
                             total_transitions, report_level);

                memcpy(&U[Ndofs * it], &xx[0], Ndofs * sizeof(int));
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
    if (t_time)
        free(t_time);
    if (sum_t_rate)
        free(sum_t_rate);
    if (t_rate)
        free(t_rate);
    if (xx)
        free(xx);
    if (update_node)
        free(update_node);

    return errcode;
}

int siminf_core(
    const int *u0, const int *irG, const int *jcG, const int *irN,
    const int *jcN, const int *prN, const double *tspan, const int tlen,
    int *U, double *data, const int *sd, const int Nn,
    const int Nc, const int Nt, const int Nobs, const int dsize,
    const int *irE, const int *jcE, const int *prE,
    const external_events *events,
    int report_level, int Nthreads, unsigned long int seed,
    const PropensityFun *t_fun, const PostTimeStepFun pts_fun,
    const ProgressFun progress, const char *strategy)
{
    int err = SIMINF_UNSUPPORTED_PARALLELIZATION;
    gsl_rng *rng = NULL;

    if (strcmp(strategy, "single") == 0) {
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);

        err = siminf_core_single(
            u0, irG, jcG, irN, jcN, prN, tspan, tlen, U, data, sd, Nn, Nc,
            Nt, Nobs, dsize, irE, jcE, prE, events, report_level, Nthreads,
            rng, t_fun, pts_fun, progress);
    }

    if (rng)
        gsl_rng_free(rng);

    return err;
}
