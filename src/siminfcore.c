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
    const int *u0, const size_t *irG, const size_t *jcG, const size_t *irN,
    const size_t *jcN, const int *prN, const double *tspan, const size_t tlen,
    int *U, double *data, const int *sd, const size_t Nn,
    const size_t Nc, const size_t Nt, const int Nobs, const size_t dsize,
    const size_t *irE, const size_t *jcE, const int *prE,
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
    size_t it = 0;
    const size_t Ndofs = Nn * Nc;

    /* Variables to handle external events */
    const int *ext_event         = events->event;
    const int *ext_time          = events->time;
    const int *ext_select        = events->select;
    const int *ext_node          = events->node;
    const int *ext_dest          = events->dest;
    const int *ext_n             = events->n;
    const double *ext_proportion = events->proportion;
    int ext_len                  = events->len;
    int ext_i                    = 0;
    ExtEventHandlerFun extfun[]  = {event_exit,
                                    event_enter,
                                    event_internal_transfer,
                                    event_external_transfer};
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
        size_t i;

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
        for (node = 0; node < Nn; node++) {
            while (t_time[node] < next_day) {
                double cum, rand, tot_rate, delta = 0.0;
                size_t i, tr = 0;

                /* a) Determine the transition that did occur (directSSA). */
                cum = t_rate[node * Nt];
                rand = gsl_rng_uniform(rng) * sum_t_rate[node];
                while (tr < Nt && rand > cum) {
                    tr++;
                    cum += t_rate[node * Nt + tr];
                }

                /* b) Update the state of the node */
                for (i = jcN[tr]; i < jcN[tr + 1]; i++) {
                    xx[node * Nc + irN[i]] += prN[i];
                    if (xx[node * Nc + irN[i]] < 0)
                        errcode = SIMINF_ERR_NEGATIVE_STATE;
                }

                /* c) Recalculate sum_t_rate[node] using dependency graph. */
                for (i = jcG[tr]; i < jcG[tr + 1]; i++) {
                    size_t j = irG[i];
                    double old = t_rate[node * Nt + j];
                    delta += (t_rate[node * Nt + j] = (*t_fun[j])(
                                  &xx[node * Nc], t_time[node],
                                  &data[node * dsize], sd[node]))
                        - old;
                }
                sum_t_rate[node] += delta;

                total_transitions++; /* counter */

                /* d) Compute time to new event for this node. */
                tot_rate = sum_t_rate[node];
                if (tot_rate > 0.0) {
                    t_time[node] = -log(1.0 - gsl_rng_uniform(rng)) / tot_rate +
                        t_time[node];
                } else {
                    t_time[node] = INFINITY;
                }

                /* e) Check for error codes. */
                if (errcode) {
                    /* Report when the error occurred. */
                    if (report_level)
                        progress(t_time[node], tspan[0], tspan[tlen - 1],
                                 total_transitions, report_level);
                    break;
                }
            }
        }

        /* Check if the exit from the while-loop was due to:
         *   a) Simulation reached the final time
         *   b) Error code. */
        if (it >= tlen || errcode)
            break;

        /* (2) Incorporate all scheduled external events. */
        while (ext_i < ext_len && tt >= ext_time[ext_i]) {
            errcode = (*extfun[ext_event[ext_i]])(irE, jcE, prE, Nc, Nobs, xx,
                                                  ext_node[ext_i], ext_dest[ext_i],
                                                  ext_select[ext_i], ext_n[ext_i],
                                                  ext_proportion[ext_i],
                                                  individuals, rng);

            /* Check for error codes. */
            if (errcode) {
                /* Report when the error occurred. */
                if (report_level)
                    progress(tt, tspan[0], tspan[tlen - 1],
                             total_transitions, report_level);
                break;
            }

            /* Indicate node and dest node for update */
            update_node[ext_node[ext_i]] = 1;
            update_node[ext_dest[ext_i]] = 1;

            ext_i++;
        }

        /* Check if the exit from the while-loop was due to error. */
        if (errcode)
            break;

        /* (3) Update the infectious pressure variable. */
        for (node = 0; node < Nn; node++) {
            if (pts_fun(&xx[node * Nc], node, tt, &data[node * dsize], sd[node]) ||
                update_node[node])
            {
                size_t i = 0;
                double delta = 0.0, old_t_rate = sum_t_rate[node];

                /* compute new transition rate only for transitions
                 * dependent on infectious pressure */
                for (; i < Nt; i++) {
                    double old = t_rate[node * Nt + i];
                    delta += (t_rate[node * Nt + i] = (*t_fun[i])(
                                  &xx[node * Nc], tt, &data[node * dsize], sd[node]))
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
    const int *u0, const size_t *irG, const size_t *jcG, const size_t *irN,
    const size_t *jcN, const int *prN, const double *tspan, const size_t tlen,
    int *U, double *data, const int *sd, const size_t Nn,
    const size_t Nc, const size_t Nt, const int Nobs, const size_t dsize,
    const size_t *irE, const size_t *jcE, const int *prE,
    const external_events *events,
    int report_level, int Nthreads, const gsl_rng *rng,
    const PropensityFun *t_fun, const PostTimeStepFun pts_fun,
    const ProgressFun progress, const char *strategy)
{
    int err = SIMINF_UNSUPPORTED_PARALLELIZATION;

    if (strcmp(strategy, "single") == 0) {
        err = siminf_core_single(
            u0, irG, jcG, irN, jcN, prN, tspan, tlen, U, data, sd, Nn, Nc,
            Nt, Nobs, dsize, irE, jcE, prE, events, report_level, Nthreads,
            rng, t_fun, pts_fun, progress);
    }

    return err;
}
