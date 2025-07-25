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

#include "SimInf.h"
#include "SimInf_internal.h"
#include "SimInf_solver.h"
#include <R_ext/Visibility.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#ifdef _OPENMP
#  include <omp.h>
#endif
#include <string.h>

/**
 * SimInf multi-model solver
 *
 * @return 0 if Ok, else error code.
 */
static int
SimInf_solver_mssm(
    SimInf_compartment_model *model,
    SimInf_scheduled_events *events,
    int Nthread)
{
    #ifdef _OPENMP
    #  pragma omp parallel num_threads(SimInf_num_threads())
    #endif
    {
        #ifdef _OPENMP
        #  pragma omp for
        #endif
        for (int i = 0; i < Nthread; i++) {
            SimInf_scheduled_events e = *&events[i];
            SimInf_compartment_model m = *&model[i];

            /* Store original pointers to the model state. */
            int *u = m.u;
            int *U = m.U;
            double *v = m.v;
            double *v_new = m.v_new;
            double *V = m.V;

            for (ptrdiff_t replicate = 0; replicate < m.Nrep && !m.error; replicate++) {
                /* Clear processed events. */
                e.events_index = 0;

                /* Move to the initial state of the model
                 * replicate. */
                m.u = &u[replicate * m.Nn * m.Nc];
                m.U = &U[replicate * m.tlen * m.Nn * m.Nc];
                m.v = &v[replicate * m.Nn * m.Nd];
                m.v_new = &v_new[replicate * m.Nn * m.Nd];
                m.V = &V[replicate * m.tlen * m.Nn * m.Nd];

                /* Initialize global time. */
                m.tt = m.tspan[0];
                m.next_unit_of_time = floor(m.tt) + 1.0;
                m.U_it = 0;
                m.V_it = 0;

                /* Initialize the transition rate for every transition
                 * and every node. Store the sum of the transition
                 * rates in each node in sum_t_rate. Moreover,
                 * initialize time in each node. */
                for (ptrdiff_t node = 0; node < m.Nn && !m.error; node++) {
                    m.sum_t_rate[node] = 0.0;
                    for (int j = 0; j < m.Nt; j++) {
                        const double rate = (*m.tr_fun[j])(
                            &m.u[node * m.Nc],
                            &m.v[node * m.Nd],
                            &m.ldata[node * m.Nld],
                            m.gdata,
                            m.tt);

                        m.t_rate[node * m.Nt + j] = rate;
                        m.sum_t_rate[node] += rate;
                        if (!R_FINITE(rate) || rate < 0.0) {
                            SimInf_print_status(
                                m.Nc,
                                &m.u[node * m.Nc],
                                m.Nd,
                                &m.v[node * m.Nd],
                                m.Nld,
                                &m.ldata[node * m.Nld],
                                (int)node,
                                m.tt,
                                rate,
                                j);
                            m.error = SIMINF_ERR_INVALID_RATE;
                        }
                    }

                    m.t_time[node] = m.tt;
                }

                /* Main loop. */
                for (;!m.error;) {
                    /* (1) Handle internal epidemiological model,
                     * continuous-time Markov chain. */
                    for (ptrdiff_t node = 0; node < m.Nn && !m.error; node++) {
                        for (;;) {
                            double cum, rand, tau, delta = 0.0;
                            int tr;

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

                            /* 1b) Determine the transition that did
                             * occur (direct SSA). */
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

                                /* No nonzero rate found, but a
                                   transition was sampled. This can
                                   happen due to floating point errors
                                   in the iterated recalculated
                                   rates. */
                                if (m.t_rate[node * m.Nt + tr] == 0.0) {
                                    /* nil event: zero out and move on */
                                    m.sum_t_rate[node] = 0.0;
                                    break;
                                }
                            }

                            /* 1c) Update the state of the node */
                            for (int j = m.jcS[tr]; j < m.jcS[tr + 1]; j++) {
                                m.u[node * m.Nc + m.irS[j]] += m.prS[j];
                                if (m.u[node * m.Nc + m.irS[j]] < 0) {
                                    SimInf_print_status(
                                        m.Nc,
                                        &m.u[node * m.Nc],
                                        m.Nd,
                                        &m.v[node * m.Nd],
                                        m.Nld,
                                        &m.ldata[node * m.Nld],
                                        (int)node,
                                        m.t_time[node],
                                        0,
                                        tr);
                                    m.error = SIMINF_ERR_NEGATIVE_STATE;
                                }
                            }

                            /* 1d) Recalculate sum_t_rate[node] using
                             * dependency graph. */
                            for (int j = m.jcG[tr]; j < m.jcG[tr + 1]; j++) {
                                const double old = m.t_rate[node * m.Nt + m.irG[j]];
                                const double rate = (*m.tr_fun[m.irG[j]])(
                                    &m.u[node * m.Nc],
                                    &m.v[node * m.Nd],
                                    &m.ldata[node * m.Nld],
                                    m.gdata,
                                    m.t_time[node]);

                                m.t_rate[node * m.Nt + m.irG[j]] = rate;
                                delta += rate - old;
                                if (!R_FINITE(rate) || rate < 0.0) {
                                    SimInf_print_status(
                                        m.Nc,
                                        &m.u[node * m.Nc],
                                        m.Nd,
                                        &m.v[node * m.Nd],
                                        m.Nld,
                                        &m.ldata[node * m.Nld],
                                        (int)node,
                                        m.t_time[node],
                                        rate,
                                        m.irG[j]);
                                    m.error = SIMINF_ERR_INVALID_RATE;
                                }
                            }
                            m.sum_t_rate[node] += delta;
                        }
                    }

                    /* (2 & 3) Incorporate all scheduled E1 and E2
                     * events */
                    SimInf_process_events(&m, &e, 1);

                    /* (4) Incorporate model specific actions after
                     * each timestep e.g. update the infectious
                     * pressure variable. Moreover, update transition
                     * rates in nodes that are indicated for update */
                    for (ptrdiff_t node = 0; node < m.Nn && !m.error; node++) {
                        const int rc = m.pts_fun(
                            &m.v_new[node * m.Nd],
                            &m.u[node * m.Nc],
                            &m.v[node * m.Nd],
                            &m.ldata[node * m.Nld],
                            m.gdata,
                            (int)node,
                            m.tt);

                        if (rc < 0) {
                            m.error = rc;
                            break;
                        } else if (rc > 0 || m.update_node[node]) {
                            /* Update transition rates */
                            double delta = 0.0;

                            for (int j = 0; j < m.Nt; j++) {
                                const double old = m.t_rate[node * m.Nt + j];
                                const double rate = (*m.tr_fun[j])(
                                    &m.u[node * m.Nc],
                                    &m.v_new[node * m.Nd],
                                    &m.ldata[node * m.Nld],
                                    m.gdata,
                                    m.tt);

                                m.t_rate[node * m.Nt + j] = rate;
                                delta += rate - old;
                                if (!R_FINITE(rate) || rate < 0.0) {
                                    SimInf_print_status(
                                        m.Nc,
                                        &m.u[node * m.Nc],
                                        m.Nd,
                                        &m.v[node * m.Nd],
                                        m.Nld,
                                        &m.ldata[node * m.Nld],
                                        (int)node,
                                        m.tt,
                                        rate,
                                        j);
                                    m.error = SIMINF_ERR_INVALID_RATE;
                                }
                            }
                            m.sum_t_rate[node] += delta;

                            m.update_node[node] = 0;
                        }
                    }

                    /* (5) The global time now equals next unit of
                     * time. */
                    m.tt = m.next_unit_of_time;
                    m.next_unit_of_time += 1.0;

                    /* (6) Store solution if tt has passed the next
                     * time in tspan. Report solution up to, but not
                     * including tt. The mssm solver always stores the
                     * solution in a dense matrix (U and/or V non-null
                     * pointers).  Copy compartment state to U */
                    while (m.U_it < m.tlen && m.tt > m.tspan[m.U_it])
                        memcpy(&m.U[(ptrdiff_t)m.Nn * (ptrdiff_t)m.Nc * m.U_it++],
                               m.u, (ptrdiff_t)m.Nn * (ptrdiff_t)m.Nc * sizeof(int));

                    /* Copy continuous state to V */
                    while (m.V_it < m.tlen && m.tt > m.tspan[m.V_it])
                        memcpy(&m.V[(ptrdiff_t)m.Nd * (ptrdiff_t)m.Ntot * m.V_it++],
                               m.v_new, (ptrdiff_t)m.Nn * (ptrdiff_t)m.Nd * sizeof(double));

                    /* Swap the pointers to the continuous state
                     * variable so that 'v' equals 'v_new'. */
                    double *v_tmp = m.v;
                    m.v = m.v_new;
                    m.v_new = v_tmp;

                    /* If the simulation has reached the final time,
                     * exit. */
                    if (m.U_it >= m.tlen)
                        break;
                }
            }

            /* Restore original pointers to the model state. */
            m.u = u;
            m.U = U;
            m.v = v;
            m.v_new = v_new;
            m.V = V;
            *&model[i] = m;
        }
    }

    /* Check for error. */
    for (int i = 0; i < Nthread; i++)
        if (model[i].error)
            return model[i].error;

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
SimInf_run_solver_mssm(
    SimInf_solver_args *args)
{
    int err = 0;
    gsl_rng *rng = NULL;
    SimInf_scheduled_events *events = NULL;
    SimInf_compartment_model *model = NULL;

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!rng) {
        err = SIMINF_ERR_ALLOC_MEMORY_BUFFER; /* #nocov */
        goto cleanup;                         /* #nocov */
    }
    gsl_rng_set(rng, args->seed);

    err = SimInf_compartment_model_create(&model, args);
    if (err)
        goto cleanup; /* #nocov */

    err = SimInf_scheduled_events_create(&events, args, rng);
    if (err)
        goto cleanup; /* #nocov */

    err = SimInf_solver_mssm(model, events, args->Nthread);

cleanup:
    gsl_rng_free(rng);
    SimInf_scheduled_events_free(events);
    SimInf_compartment_model_free(model);

    return err;
}
