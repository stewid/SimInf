/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015 Pavol Bauer
 *  Copyright (C) 2017 - 2018 Robin Eriksson
 *  Copyright (C) 2015 - 2018 Stefan Engblom
 *  Copyright (C) 2015 - 2018 Stefan Widgren
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "SimInf.h"
#include "SimInf_solver_ssa.h"

/**
 * Siminf solver
 *
 * @return 0 if Ok, else error code.
 */
static int SimInf_solver_ssa(
    SimInf_compartment_model *model, SimInf_scheduled_events *events,
    int *uu, int *update_node, int Nthread)
{
    int k;

    #pragma omp parallel
    {
        int i;

        #pragma omp for
        for (i = 0; i < Nthread; i++) {
            int node;
            SimInf_compartment_model m = *&model[i];

            /* Initialize the transition rate for every transition and
             * every node. Store the sum of the transition rates in
             * each node in sum_t_rate. Moreover, initialize time in
             * each node. */
            for (node = 0; node < m.Nn; node++) {
                int j;

                m.sum_t_rate[node] = 0.0;
                for (j = 0; j < m.Nt; j++) {
                    const double rate = (*m.tr_fun[j])(
                            &m.u[node * m.Nc], &m.v[node * m.Nd],
                            &m.ldata[node * m.Nld], m.gdata, m.tt);

                    m.t_rate[node * m.Nt + j] = rate;
                    m.sum_t_rate[node] += rate;
                    if (!isfinite(rate) || rate < 0.0)
                        m.errcode = SIMINF_ERR_INVALID_RATE;
                }

                m.t_time[node] = m.tt;
            }

            *&model[i] = m;
        }
    }

    /* Check for error during initialization. */
    for (k = 0; k < Nthread; k++)
        if (model[k].errcode)
            return model[k].errcode;

    /* Main loop. */
    for (;;) {
        #pragma omp parallel
        {
            int i;

            #pragma omp for
            for (i = 0; i < Nthread; i++) {
                int node;
                SimInf_scheduled_events e = *&events[i];
                SimInf_compartment_model m = *&model[i];

                /* (1) Handle internal epidemiological model,
                 * continuous-time Markov chain. */
                for (node = 0; node < m.Nn && !m.errcode; node++) {
                    for (;;) {
                        double cum, rand, tau, delta = 0.0;
                        int j, tr;

                        /* 1a) Compute time to next event for this
                         * node. */
                        if (m.sum_t_rate[node] <= 0.0) {
                            m.t_time[node] = m.next_day;
                            break;
                        }
                        tau = -log(gsl_rng_uniform_pos(e.rng)) /
                            m.sum_t_rate[node];
                        if ((tau + m.t_time[node]) >= m.next_day) {
                            m.t_time[node] = m.next_day;
                            break;
                        }
                        m.t_time[node] += tau;

                        /* 1b) Determine the transition that did occur
                         * (direct SSA). */
                        rand = gsl_rng_uniform_pos(e.rng) * m.sum_t_rate[node];
                        for (tr = 0, cum = m.t_rate[node * m.Nt];
                             tr < m.Nt && rand > cum;
                             tr++, cum += m.t_rate[node * m.Nt + tr]);

                        /* Elaborate floating point fix: */
                        if (tr >= m.Nt)
                            tr = m.Nt - 1;
                        if (m.t_rate[node * m.Nt + tr] == 0.0) {
                            /* Go backwards and try to find first
                             * nonzero transition rate */
                            for ( ; tr > 0 && m.t_rate[node * m.Nt + tr] == 0.0; tr--);

                            /* No nonzero rate found, but a transition
                               was sampled. This can happen due to
                               floating point errors in the iterated
                               recalculated rates. */
                            if (m.t_rate[node * m.Nt + tr] == 0.0) {
                                /* nil event: zero out and move on */
                                m.sum_t_rate[node] = 0.0;
                                break;
                            }
                        }

                        /* 1c) Update the state of the node */
                        for (j = m.jcS[tr]; j < m.jcS[tr + 1]; j++) {
                            m.u[node * m.Nc + m.irS[j]] += m.prS[j];
                            if (m.u[node * m.Nc + m.irS[j]] < 0)
                                m.errcode = SIMINF_ERR_NEGATIVE_STATE;
                        }

                        /* 1d) Recalculate sum_t_rate[node] using
                         * dependency graph. */
                        for (j = m.jcG[tr]; j < m.jcG[tr + 1]; j++) {
                            const double old = m.t_rate[node * m.Nt + m.irG[j]];
                            const double rate = (*m.tr_fun[m.irG[j]])(
                                &m.u[node * m.Nc], &m.v[node * m.Nd],
                                &m.ldata[node * m.Nld], m.gdata,
                                m.t_time[node]);

                            m.t_rate[node * m.Nt + m.irG[j]] = rate;
                            delta += rate - old;
                            if (!isfinite(rate) || rate < 0.0)
                                m.errcode = SIMINF_ERR_INVALID_RATE;
                        }
                        m.sum_t_rate[node] += delta;
                    }
                }

                *&events[i] = e;
                *&model[i] = m;

                /* (2) Incorporate all scheduled E1 events */
                SimInf_process_E1_events(&model[i], &events[i], uu, update_node);
            }

            #pragma omp barrier

            #pragma omp master
            {
                /* (3) Incorporate all scheduled E2 events */
                SimInf_process_E2_events(model, events, uu, update_node);
            }

            #pragma omp barrier

            #pragma omp for
            for (i = 0; i < Nthread; i++) {
                int node;
                SimInf_compartment_model m = *&model[i];

                /* (4) Incorporate model specific actions after each
                 * timestep e.g. update the infectious pressure
                 * variable. Moreover, update transition rates in
                 * nodes that are indicated for update */
                for (node = 0; node < m.Nn; node++) {
                    const int rc = m.pts_fun(
                        &m.v_new[node * m.Nd], &m.u[node * m.Nc],
                        &m.v[node * m.Nd], &m.ldata[node * m.Nld],
                        m.gdata, m.Ni + node, m.tt);

                    if (rc < 0) {
                        m.errcode = rc;
                        break;
                    } else if (rc > 0 || m.update_node[node]) {
                        /* Update transition rates */
                        int j = 0;
                        double delta = 0.0;

                        for (; j < m.Nt; j++) {
                            const double old = m.t_rate[node * m.Nt + j];
                            const double rate = (*m.tr_fun[j])(
                                &m.u[node * m.Nc], &m.v_new[node * m.Nd],
                                &m.ldata[node * m.Nld], m.gdata, m.tt);

                            m.t_rate[node * m.Nt + j] = rate;
                            delta += rate - old;
                            if (!isfinite(rate) || rate < 0.0)
                                m.errcode = SIMINF_ERR_INVALID_RATE;
                        }
                        m.sum_t_rate[node] += delta;

                        m.update_node[node] = 0;
                    }
                }

                /* (5) The global time now equals next_day. */
                m.tt = m.next_day;
                m.next_day += 1.0;

                /* (6) Store solution if tt has passed the next time
                 * in tspan. Report solution up to, but not including
                 * tt. The default is to store the solution in a dense
                 * matrix (U and/or V non-null pointers) (6a).
                 * However, it is possible to store the solution in a
                 * sparse matrix. In that case, the solution is stored
                 * outside the 'pragma omp parallel' statement (6b). */
                /* 6a) Handle the case where the solution is stored in
                 * a dense matrix */
                /* Copy compartment state to U */
                while (m.U && m.U_it < m.tlen && m.tt > m.tspan[m.U_it])
                    memcpy(&m.U[m.Nc * ((m.Ntot * m.U_it++) + m.Ni)],
                           m.u, m.Nn * m.Nc * sizeof(int));
                /* Copy continuous state to V */
                while (m.V && m.V_it < m.tlen && m.tt > m.tspan[m.V_it])
                    memcpy(&m.V[m.Nd * ((m.Ntot * m.V_it++) + m.Ni)],
                           m.v_new, m.Nn * m.Nd * sizeof(double));

                *&model[i] = m;
            }
        }

        /* 6b) Handle the case where the solution is stored in a sparse
         * matrix */
        SimInf_store_solution_sparse(model);

        /* Swap the pointers to the continuous state variable so that
         * 'v' equals 'v_new'. Moreover, check for error. */
        for (k = 0; k < Nthread; k++) {
            double *v_tmp = model[k].v;
            model[k].v = model[k].v_new;
            model[k].v_new = v_tmp;
            if (model[k].errcode)
                return model[k].errcode;
        }

        /* If the simulation has reached the final time, exit. */
        if (model[0].U_it >= model[0].tlen)
            break;
    }

    return 0;
}

/**
 * Initialize and run siminf solver
 *
 * @param args Structure with data for the solver.
 * @return 0 if Ok, else error code.
 */
int SimInf_run_solver_ssa(SimInf_solver_args *args)
{
    int error = 0, i;
    gsl_rng *rng = NULL;
    SimInf_scheduled_events *events = NULL;
    SimInf_compartment_model *model = NULL;
    int *uu = NULL, *update_node = NULL;
    double *vv_1 = NULL, *vv_2 = NULL;

    /* Set compartment state to the initial state. */
    uu = malloc(args->Nn * args->Nc * sizeof(int));
    if (!uu) {
        error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    memcpy(uu, args->u0, args->Nn * args->Nc * sizeof(int));

    /* Copy u0 to either U[, 1] or U_sparse[, 1] */
    if (args->U) {
        memcpy(args->U, args->u0, args->Nn * args->Nc * sizeof(int));
    } else {
        for (i = args->jcU[0]; i < args->jcU[1]; i++)
            args->prU[i] = args->u0[args->irU[i]];
    }

    /* Set continuous state to the initial state in each node. */
    vv_1 = malloc(args->Nn * args->Nd * sizeof(double));
    vv_2 = malloc(args->Nn * args->Nd * sizeof(double));
    if (!vv_1 || !vv_2) {
        error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    memcpy(vv_1, args->v0, args->Nn * args->Nd * sizeof(double));

    /* Copy v0 to either V[, 1] or V_sparse[, 1] */
    if (args->V) {
        memcpy(args->V, args->v0, args->Nn * args->Nd * sizeof(double));
    } else {
        for (i = args->jcV[0]; i < args->jcV[1]; i++)
            args->prV[i] = args->v0[args->irV[i]];
    }

    /* Setup vector to keep track of nodes that must be updated due to
     * scheduled events */
    update_node = calloc(args->Nn, sizeof(int));
    if (!update_node) {
        error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!rng) {
        error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    gsl_rng_set(rng, args->seed);

    error = SimInf_compartment_model_create(
        &model, args, uu, vv_1, vv_2, update_node);
    if (error)
        goto cleanup;

    error = SimInf_scheduled_events_create(&events, args, rng);
    if (error)
        goto cleanup;

    error = SimInf_solver_ssa(model, events, uu, update_node, args->Nthread);

cleanup:
    if (uu)
        free(uu);

    if (vv_1)
        free(vv_1);

    if (vv_2)
        free(vv_2);

    if (update_node)
        free(update_node);

    if (rng)
        gsl_rng_free(rng);

    SimInf_scheduled_events_free(events, args->Nthread);
    SimInf_compartment_model_free(model, args->Nthread);

    return error;
}
