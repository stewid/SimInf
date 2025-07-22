/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
 * Copyright (C) 2015 -- 2025 Stefan Widgren
 *
 * SimInf is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SimInf is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <R_ext/Visibility.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "SimInf.h"
#include "misc/SimInf_openmp.h"
#include "SimInf_solver.h"

/**
 * Siminf solver
 *
 * @return 0 if Ok, else error code.
 */
static int
SimInf_solver_ssm(
    SimInf_compartment_model *model,
    SimInf_scheduled_events *events)
{
    int Nthread = model->Nthread;
    int k;

    #ifdef _OPENMP
    #  pragma omp parallel num_threads(SimInf_num_threads())
    #endif
    {
        int i;

        #ifdef _OPENMP
        #  pragma omp for
        #endif
        for (i = 0; i < Nthread; i++) {
            SimInf_compartment_model m = *&model[i];

            /* Initialize the transition rate for every transition and
             * every node. Store the sum of the transition rates in
             * each node in sum_t_rate. Moreover, initialize time in
             * each node. */
            for (int node = 0; node < m.Nn; node++) {
                m.sum_t_rate[node] = 0.0;
                for (int j = 0; j < m.Nt; j++) {
                    const double rate = (*m.tr_fun[j])(
                            &m.u[node * m.Nc], &m.v[node * m.Nd],
                            &m.ldata[node * m.Nld], m.gdata, m.tt);

                    m.t_rate[node * m.Nt + j] = rate;
                    m.sum_t_rate[node] += rate;
                    if (!R_FINITE(rate) || rate < 0.0) {
                        SimInf_print_status(m.Nc, &m.u[node * m.Nc],
                                            m.Nd, &m.v[node * m.Nd],
                                            m.Nld, &m.ldata[node * m.Nld],
                                            m.Ni + node, m.tt, rate, j);
                        m.error = SIMINF_ERR_INVALID_RATE;
                    }
                }

                m.t_time[node] = m.tt;
            }

            *&model[i] = m;
        }
    }

    /* Check for error during initialization. */
    for (k = 0; k < Nthread; k++)
        if (model[k].error)
            return model[k].error;

    /* Main loop. */
    for (;;) {
        #ifdef _OPENMP
        #  pragma omp parallel num_threads(SimInf_num_threads())
        #endif
        {
            int i;

            #ifdef _OPENMP
            #  pragma omp for
            #endif
            for (i = 0; i < Nthread; i++) {
                SimInf_scheduled_events e = *&events[i];
                SimInf_compartment_model m = *&model[i];

                /* (1) Handle internal epidemiological model,
                 * continuous-time Markov chain. */
                for (int node = 0; node < m.Nn && !m.error; node++) {
                    for (;;) {
                        double cum, rand, tau, delta = 0.0;
                        int j, tr;

                        /* 1a) Compute time to next event for this
                         * node. */
                        if (m.sum_t_rate[node] <= 0.0) {
                            m.t_time[node] = m.next_unit_of_time;
                            break;
                        }
                        tau = -log(gsl_rng_uniform_pos(e.rng)) /
                            m.sum_t_rate[node];
                        if ((tau + m.t_time[node]) >= m.next_unit_of_time) {
                            m.t_time[node] = m.next_unit_of_time;
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
                            if (m.u[node * m.Nc + m.irS[j]] < 0) {
                                SimInf_print_status(m.Nc, &m.u[node * m.Nc],
                                                    m.Nd, &m.v[node * m.Nd],
                                                    m.Nld, &m.ldata[node * m.Nld],
                                                    m.Ni + node, m.t_time[node],
                                                    0, tr);
                                m.error = SIMINF_ERR_NEGATIVE_STATE;
                            }
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
                            if (!R_FINITE(rate) || rate < 0.0) {
                                SimInf_print_status(m.Nc, &m.u[node * m.Nc],
                                                    m.Nd, &m.v[node * m.Nd],
                                                    m.Nld, &m.ldata[node * m.Nld],
                                                    m.Ni + node, m.t_time[node],
                                                    rate, m.irG[j]);
                                m.error = SIMINF_ERR_INVALID_RATE;
                            }
                        }
                        m.sum_t_rate[node] += delta;
                    }
                }

                *&events[i] = e;
                *&model[i] = m;

                /* (2) Incorporate all scheduled E1 events */
                SimInf_process_events(&model[i], &events[i], 0);
            }

            #ifdef _OPENMP
            #  pragma omp barrier
            #endif

            #ifdef _OPENMP
            #  pragma omp master
            #endif
            {
                /* (3) Incorporate all scheduled E2 events */
                SimInf_process_events(model, events, 1);
            }

            #ifdef _OPENMP
            #  pragma omp barrier
            #endif

            #ifdef _OPENMP
            #  pragma omp for
            #endif
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
                        m.error = rc;
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
                            if (!R_FINITE(rate) || rate < 0.0) {
                                SimInf_print_status(m.Nc, &m.u[node * m.Nc],
                                                    m.Nd, &m.v[node * m.Nd],
                                                    m.Nld, &m.ldata[node * m.Nld],
                                                    m.Ni + node, m.tt, rate, j);
                                m.error = SIMINF_ERR_INVALID_RATE;
                            }
                        }
                        m.sum_t_rate[node] += delta;

                        m.update_node[node] = 0;
                    }
                }

                /* (5) The global time now equals next unit of time. */
                m.tt = m.next_unit_of_time;
                m.next_unit_of_time += 1.0;

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
            if (model[k].error)
                return model[k].error;
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
attribute_hidden
int
SimInf_run_solver_ssm(
    SimInf_solver_args *args)
{
    int error = 0;
    gsl_rng *rng = NULL;
    SimInf_scheduled_events *events = NULL;
    SimInf_compartment_model *model = NULL;

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!rng) {
        error = SIMINF_ERR_ALLOC_MEMORY_BUFFER; /* #nocov */
        goto cleanup;                           /* #nocov */
    }
    gsl_rng_set(rng, args->seed);

    error = SimInf_compartment_model_create(&model, args);
    if (error)
        goto cleanup; /* #nocov */

    error = SimInf_scheduled_events_create(&events, args, rng);
    if (error)
        goto cleanup; /* #nocov */

    error = SimInf_solver_ssm(model, events);

cleanup:
    gsl_rng_free(rng);
    SimInf_scheduled_events_free(events);
    SimInf_compartment_model_free(model);

    return error;
}
