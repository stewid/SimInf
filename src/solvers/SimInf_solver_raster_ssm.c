/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
 * Copyright (C) 2015 -- 2026 Stefan Widgren
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
#include <R_ext/Visibility.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#ifdef _OPENMP
#  include <omp.h>
#endif
#include <string.h>

typedef struct SimInf_raster_model
{
    /*** Data vectors for propensities ***/
    double *sum_cell_rate; /**< Vector of length Ncells with the sum
                            *   of propensities in every cell. */
    double *cell_rate;     /**< Transition rate matrix (Nt X Ncells)
                            *   with all propensities for state
                            *   transitions in cells. */

    /*** Binary (min)heap. ***/
    int *cells;
    int *heap;

    /*** Keep track of time ***/
    double *cell_time;   /**< Time for next event in each cell. */

    /*** Callbacks ***/
    TRRasterFun *tr_fun;  /**< Vector of function pointers to
                           *   transition rate functions. */

    /*** Data vectors for the cells ***/
    const int *raster;   /**< The raster data vector where raster[i]
                          *   gives the data value for cell #i. */
    int *cell_u;         /**< Vector with the count of each compartment in
                          *   each cell. */

    /*** Constants ***/
    int nrow;       /**< Number of rows. */
    int ncol;       /**< Number of cols. */
    int cell_Nc;    /**< Number of compartments in each cell. */
    int cell_i;     /**< Index to the cell compartment in each
                     *   node. */
} SimInf_raster_model;

/**
 * Free allocated memory for a raster model.
 *
 * @param model the data structure to free.
 */
static void
SimInf_raster_model_free(
    SimInf_raster_model *model)
{
    if (model) {
        /* Free data vectors for propensities. */
        free(model->sum_cell_rate);
        model->sum_cell_rate = NULL;
        free(model->cell_rate);
        model->cell_rate = NULL;

        /* Free data for binary (min)heap. */
        free(model->cells);
        model->cells = NULL;
        free(model->heap);
        model->heap = NULL;

        /* Free data to keep track of time. */
        free(model->cell_time);
        model->cell_time = NULL;

        /* Free data vector to keep track of the states in the
         * cells. */
        free(model->cell_u);
        model->cell_u = NULL;
    }
}

/**
 * Create and initialize data to process a raster model. The generated
 * data structure must be freed by the user.
 *
 * @param out the resulting data structure.
 * @param args structure with data for the solver.
 * @param rng random number generator
 * @return 0 or an error code
 */
static int
SimInf_raster_model_create(
    SimInf_raster_model **out,
    SimInf_solver_args const *args)
{
    SimInf_raster_model *model = NULL;

    model = calloc(1, sizeof(SimInf_raster_model));
    if (!model)
        goto on_error;

    /* Callbacks */
    model->tr_fun = args->tr_raster_fun;

    /* Constants. */
    model->nrow = args->nrow;
    model->ncol = args->ncol;
    if (model->nrow < 1 || model->ncol < 1)
        goto on_error;

    /* Number of compartments in each cell. */
    model->cell_Nc = args->cell_Nc;
    if (model->cell_Nc < 0)
        goto on_error;

    /* Index to the cell compartment in each node. */
    model->cell_i = args->cell_i;
    if (model->cell_i < 0 || model->cell_i >= args->Nc)
        goto on_error;

    /* Raster data. */
    model->raster = args->raster;

    /* Number of compartments in each cell. */
    model->cell_Nc = args->cell_Nc;
    if (model->cell_Nc < 0)
        goto on_error;

    /* Allocate memory for the cell compartment state and set it to
     * zero, the initial state. */
    model->cell_u = calloc(model->nrow * model->ncol * model->cell_Nc,
                           sizeof(int));
    if (!model->cell_u)
        goto on_error;

    /* Binary (min)heap. */
    model->cells = malloc(model->nrow * model->ncol * sizeof(int));
    if (!model->cells)
        goto on_error;

    model->heap = malloc(model->nrow * model->ncol * sizeof(int));
    if (!model->heap)
        goto on_error;

    model->cell_time = malloc(model->nrow * model->ncol * sizeof(double));
    if (!model->cell_time)
        goto on_error;

    /* Create transition rate matrix (Nt X Ncells) and total rate
     * vector. In cell_rate we store all propensities for state
     * transitions, and in sum_cell_rate the sum of propensities in
     * every cell. */
    model->cell_rate = calloc(args->Nt * model->nrow * model->ncol,
                              sizeof(double));
    if (!model->cell_rate)
        goto on_error;
    model->sum_cell_rate = calloc(model->nrow * model->ncol,
                                  sizeof(double));
    if (!model->sum_cell_rate)
        goto on_error;

    *out = model;
    return 0;

on_error:
    SimInf_raster_model_free(model);
    return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
}

/**
 * Siminf solver
 *
 * @return 0 if Ok, else error code.
 */
static int
SimInf_solver_raster_ssm(
    SimInf_compartment_model *model,
    SimInf_scheduled_events *events,
    SimInf_raster_model *raster)
{
    bool done = false;
    int Nthread = model->Nthread;
    int Ncells = raster->nrow * raster->ncol;

#ifdef _OPENMP
#  pragma omp parallel num_threads(SimInf_num_threads())
#endif
    {
#ifdef _OPENMP
#  pragma omp for schedule(static)
#endif
        for (int i = 0; i < Nthread; i++) {
            SimInf_compartment_model m = *&model[i];

            /* Initialize the transition rate for every transition and
             * every node. Store the sum of the transition rates in
             * each node in sum_t_rate. Moreover, initialize time in
             * each node. */
            for (ptrdiff_t node = 0; node < m.Nn; node++) {
                /* Note that the cell location in R is one-based,
                 * therefore, decrement the cell with one. */
                const int cell = m.u[node * m.Nc + raster->cell_i] - 1;
                const int *cell_u = &raster->cell_u[cell * raster->cell_Nc];

                m.sum_t_rate[node] = 0.0;
                for (int j = 0; j < m.Nt; j++) {
                    const double rate = (*raster->tr_fun[j]) (
                        raster->raster,
                        raster->nrow,
                        raster->ncol,
                        cell_u,
                        &m.u[node * m.Nc],
                        &m.v[node * m.Nd],
                        &m.ldata[node * m.Nld],
                        m.gdata,
                        m.tt);

                    m.t_rate[node * m.Nt + j] = rate;
                    m.sum_t_rate[node] += rate;
                    if (!R_FINITE(rate) || rate < 0.0) {
                        SimInf_print_status(m.Nc, &m.u[node * m.Nc],
                                            m.Nd, &m.v[node * m.Nd],
                                            m.Nld, &m.ldata[node * m.Nld],
                                            (int) (m.Ni + node), m.tt, rate, j);
                        m.error = SIMINF_ERR_INVALID_RATE;
                    }
                }

                m.t_time[node] = m.tt;
            }

            *&model[i] = m;
        }

#ifdef _OPENMP
#  pragma omp single
#endif
        {
            /* Check for error during initialization. */
            for (int i = 0; i < Nthread; i++) {
                if (model[i].error)
                    done = true;
            }

            /* Next, initialize the time to the next event in every
             * cell. Initially, all transition rates are zero. */
            if (!done) {
                for (int cell = 0; cell < Ncells; cell++) {
                    raster->cell_time[cell] = R_PosInf;
                    raster->heap[cell] = cell;
                    raster->cells[cell] = cell;
                }

                initialize_heap(
                    raster->cell_time,
                    raster->cells,
                    raster->heap,
                    Ncells);
            }
        }

        /* Main loop. */
        while (!done) {
#ifdef _OPENMP
#  pragma omp for schedule(static)
#endif
            for (int i = 0; i < Nthread; i++) {
                SimInf_scheduled_events e = *&events[i];
                SimInf_compartment_model m = *&model[i];

                /* Handle the continuous-time Markov chain for the
                 * epidemiological compartment model. */
                for (ptrdiff_t node = 0; node < m.Nn && !m.error; node++) {
                    while (true) {
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

                        /* Determine the transition that did occur
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
                            for (;
                                 tr > 0
                                 && m.t_rate[node * m.Nt + tr] == 0.0; tr--);

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

                        /* Update the state of the node */
                        for (int j = m.jcS[tr]; j < m.jcS[tr + 1]; j++) {
                            m.u[node * m.Nc + m.irS[j]] += m.prS[j];
                            if (m.u[node * m.Nc + m.irS[j]] < 0) {
                                SimInf_print_status(m.Nc,
                                                    &m.u[node * m.Nc],
                                                    m.Nd,
                                                    &m.v[node * m.Nd],
                                                    m.Nld,
                                                    &m.ldata[node * m.Nld],
                                                    (int) (m.Ni + node),
                                                    m.t_time[node], 0, tr);
                                m.error = SIMINF_ERR_NEGATIVE_STATE;
                            }
                        }

                        /* Recalculate sum_t_rate[node] using the
                         * dependency graph.  Note that the cell
                         * location in R is one-based, therefore,
                         * decrement the cell with one. */
                        const int cell = m.u[node * m.Nc + raster->cell_i] - 1;
                        const int *cell_u = &raster->cell_u[cell * raster->cell_Nc];
                        for (int j = m.jcG[tr]; j < m.jcG[tr + 1]; j++) {
                            const double old = m.t_rate[node * m.Nt + m.irG[j]];
                            const double rate =
                                (*raster->tr_fun[m.irG[j]]) (
                                    raster->raster,
                                    raster->nrow,
                                    raster->ncol,
                                    cell_u,
                                    &m.u[node * m.Nc],
                                    &m.v[node * m.Nd],
                                    &m.ldata[node * m.Nld],
                                    m.gdata,
                                    m.t_time[node]);

                            m.t_rate[node * m.Nt + m.irG[j]] = rate;
                            delta += rate - old;
                            if (!R_FINITE(rate) || rate < 0.0) {
                                SimInf_print_status(m.Nc,
                                                    &m.u[node * m.Nc],
                                                    m.Nd,
                                                    &m.v[node * m.Nd],
                                                    m.Nld,
                                                    &m.ldata[node * m.Nld],
                                                    (int) (m.Ni + node),
                                                    m.t_time[node], rate,
                                                    m.irG[j]);
                                m.error = SIMINF_ERR_INVALID_RATE;
                            }
                        }
                        m.sum_t_rate[node] += delta;
                    }
                }

                *&events[i] = e;
                *&model[i] = m;

                /* Incorporate all scheduled E1 events, i.e., those
                 * events that operate on the compartments of a single
                 * node E1 = {enter, internal transfer, exit}. */
                SimInf_process_events(&model[i], &events[i], 0);
            }

#ifdef _OPENMP
#  pragma omp single
#endif
            {
                /* Incorporate all scheduled E2 events, i.e., those
                 * events that operate on the compartments of two
                 * nodes E2 = {external transfer}. */
                SimInf_process_events(model, events, 1);
            }

#ifdef _OPENMP
#  pragma omp for schedule(static)
#endif
            for (int i = 0; i < Nthread; i++) {
                SimInf_compartment_model m = *&model[i];

                /* Incorporate model specific actions after each
                 * timestep e.g. update the infectious pressure
                 * variable. Moreover, update transition rates in
                 * nodes that are indicated for update */
                for (ptrdiff_t node = 0; node < m.Nn; node++) {
                    const int rc =
                        m.pts_fun(&m.v_new[node * m.Nd], &m.u[node * m.Nc],
                                  &m.v[node * m.Nd],
                                  &m.ldata[node * m.Nld],
                                  m.gdata, (int) (m.Ni + node), m.tt);

                    if (rc < 0) {
                        m.error = rc;
                        break;
                    } else if (rc > 0 || m.update_node[node]) {
                        /* Update the transition rates.  Note that the
                         *  cell location in R is one-based,
                         *  therefore, decrement the cell with one. */
                        const int cell = m.u[node * m.Nc + raster->cell_i] - 1;
                        const int *cell_u = &raster->cell_u[cell * raster->cell_Nc];
                        double delta = 0.0;

                        for (int j = 0; j < m.Nt; j++) {
                            const double old = m.t_rate[node * m.Nt + j];
                            const double rate =
                                (*raster->tr_fun[j]) (
                                    raster->raster,
                                    raster->nrow,
                                    raster->ncol,
                                    cell_u,
                                    &m.u[node * m.Nc],
                                    &m.v_new[node * m.Nd],
                                    &m.ldata[node * m.Nld],
                                    m.gdata,
                                    m.tt);

                            m.t_rate[node * m.Nt + j] = rate;
                            delta += rate - old;
                            if (!R_FINITE(rate) || rate < 0.0) {
                                SimInf_print_status(m.Nc,
                                                    &m.u[node * m.Nc],
                                                    m.Nd,
                                                    &m.v[node * m.Nd],
                                                    m.Nld,
                                                    &m.ldata[node * m.Nld],
                                                    (int) (m.Ni + node),
                                                    m.tt, rate, j);
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
                while (m.U && m.U_it < m.tlen && m.tt > m.tspan[m.U_it]) {
                    memcpy(&m.U[(ptrdiff_t) m.Nc *
                                ((m.Ntot * m.U_it++) + m.Ni)], m.u,
                           (ptrdiff_t) m.Nn * (ptrdiff_t) m.Nc * sizeof(int));
                }

                /* Copy continuous state to V */
                while (m.V && m.V_it < m.tlen && m.tt > m.tspan[m.V_it]) {
                    const ptrdiff_t v_len = (ptrdiff_t) m.Nn * (ptrdiff_t) m.Nd *
                        sizeof(double);
                    if (v_len > 0) {
                        memcpy(&m.V[(ptrdiff_t) m.Nd * ((m.Ntot * m.V_it) + m.Ni)],
                               m.v_new,
                               v_len);
                    }
                    m.V_it++;
                }

                *&model[i] = m;
            }

#ifdef _OPENMP
#  pragma omp single
#endif
            {
                /* 6b) Handle the case where the solution is stored in
                 * a sparse matrix */
                SimInf_store_solution_sparse(model);

                /* Swap the pointers to the continuous state variable
                 * so that 'v' equals 'v_new'. Moreover, check for
                 * error. */
                for (int i = 0; i < Nthread; i++) {
                    double *v_tmp = model[i].v;
                    model[i].v = model[i].v_new;
                    model[i].v_new = v_tmp;
                    if (model[i].error)
                        done = true;
                }

                /* If the simulation has reached the final time,
                 * exit. */
                if (model[0].U_it >= model[0].tlen)
                    done = true;
            }
        } /* End of while(!done) */
    }  /* end parallel region */

    /* Check if there is any error during the simulation. */
    for (int i = 0; i < Nthread; i++) {
        if (model[i].error)
            return model[i].error;
    }

    return 0;
}

/**
 * Initialize and run siminf solver
 *
 * @param args Structure with data for the solver.
 * @return 0 if Ok, else error code.
 */
attribute_hidden int
SimInf_run_solver_raster(
    SimInf_solver_args *args)
{
    int err = 0;
    gsl_rng *rng = NULL;
    SimInf_scheduled_events *events = NULL;
    SimInf_compartment_model *model = NULL;
    SimInf_raster_model* raster = NULL;

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!rng) {
        err = SIMINF_ERR_ALLOC_MEMORY_BUFFER;   /* #nocov */
        goto cleanup;           /* #nocov */
    }
    gsl_rng_set(rng, args->seed);

    err = SimInf_compartment_model_create(&model, args);
    if (err)
        goto cleanup;           /* #nocov */

    err = SimInf_scheduled_events_create(&events, args, rng);
    if (err)
        goto cleanup;           /* #nocov */

    err = SimInf_raster_model_create(&raster, args);
    if (err)
        goto cleanup;           /* #nocov */

    err = SimInf_solver_raster_ssm(model, events, raster);

cleanup:
    gsl_rng_free(rng);
    SimInf_scheduled_events_free(events);
    SimInf_compartment_model_free(model);
    SimInf_raster_model_free(raster);

    return err;
}
