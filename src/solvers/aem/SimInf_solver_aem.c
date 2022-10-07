/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
 * Copyright (C) 2015 -- 2022 Stefan Widgren
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
#include <float.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "SimInf.h"
#include "misc/SimInf_openmp.h"
#include "SimInf_solver_aem.h"
#include "misc/binheap.h"

/**
 * Structure to hold AEM solver specific data/arguments for simulation.
 */
typedef struct SimInf_aem_arguments
{
    gsl_rng **rng_vec;   /**< The random number generator. */

    int *reactHeap;      /**< Binary heap storing all reaction events */
    int *reactNode;
    double *reactTimes;
    double *reactInf;
    int reactHeapSize;

} SimInf_aem_arguments;

/**
 * Calculate update of waiting times, including sleeping times.
 *
 * @param time, time vector.
 * @param infTime, vector of with information about time before rate was set to zero.
 * @param tt, time until next time step.
 * @param old_rate, rate before recalibration.
 * @param new_rate, rate after event.
 * @param rng, current value of rng in heap.
 */
static void
calcTimes(double* time, double* infTime, double tt, double old_rate,
          double new_rate, gsl_rng* rng)
{
    double oldtime = time[0];

    if (isinf(oldtime)) {
        if (infTime[0] == 0.0) // Waking up first time
            time[0] = -log(gsl_rng_uniform_pos(rng)) / new_rate + tt;
        else if (new_rate > 0.0)  // Waking up the 2nd..nth time
            time[0] = tt + (infTime[0] / new_rate);
    } else if (new_rate >= DBL_MIN) {
        if (oldtime == tt) // Regular update of current event
            time[0] = -log(gsl_rng_uniform_pos(rng)) / new_rate + tt;
        else  // Regular update of dependent events (rescaling)
            time[0] = ((old_rate / new_rate) * (oldtime - tt)) + tt;
    } else { // Next event time set to infinity
        infTime[0] = (oldtime - tt) * old_rate;
        time[0] = INFINITY;
    }
}

/**
 * Siminf solver
 *
 * @return 0 if Ok, else error code.
 */
static int SimInf_solver_aem(
    SimInf_compartment_model *model,
    SimInf_aem_arguments *method,
    SimInf_scheduled_events *events,
    int Nthread)
{
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
            int node;
            SimInf_compartment_model sa = *&model[i];
            SimInf_aem_arguments ma = *&method[i];

            /* Initialize the transition rate for every transition and
             * every node. */

	    /* Calculate the propensity for every reaction*/
	    for (node = 0; node < sa.Nn; node++) {
                int j;
                for (j = 0; j < sa.Nt; j++){
                    const double rate = (*sa.tr_fun[j])(&sa.u[node * sa.Nc],
                                                        &sa.v[node * sa.Nd],
                                                        &sa.ldata[node * sa.Nld],
                                                        sa.gdata,
                                                        sa.tt);
                    sa.t_rate[node * sa.Nt + j] = rate;

                    if (!R_FINITE(rate) || rate < 0.0) {
                        SimInf_print_status(sa.Nc, &sa.u[node * sa.Nc],
                                            sa.Nd, &sa.v[node * sa.Nd],
                                            sa.Nld, &sa.ldata[node * sa.Nld],
                                            sa.Ni + node, sa.tt, rate, j);
                        sa.error = SIMINF_ERR_INVALID_RATE;
                    }

                    /* calculate time until next transition j event */
                    ma.reactTimes[sa.Nt*node+j] =  -log(gsl_rng_uniform_pos(ma.rng_vec[sa.Nt*node+j]))/rate + sa.tt;
                    if (ma.reactTimes[sa.Nt*node+j] <= 0.0)
                        ma.reactTimes[sa.Nt*node+j] = INFINITY;

                    ma.reactHeap[sa.Nt*node+j] = ma.reactNode[sa.Nt*node+j] = j;
                }

                /* Initialize reaction heap */
                initialize_heap(&ma.reactTimes[sa.Nt*node], &ma.reactNode[sa.Nt*node],
                                &ma.reactHeap[sa.Nt*node], ma.reactHeapSize);
                sa.t_time[node] = sa.tt;
	    }

	    *&model[i] = sa;
            *&method[i] = ma;
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
                int node;
                SimInf_compartment_model sa = *&model[i];
                SimInf_aem_arguments ma = *&method[i];

                /* (1) Handle internal epidemiological model,
                 * continuous-time Markov chain. */
                for (node = 0; node < sa.Nn && !sa.error; node++) {
                    for (;;) {
                        int ii,j,tr;
                        double old_t_rate,rate;

                        /* 1a) Step time forward until next event */
                        sa.t_time[node] = ma.reactTimes[sa.Nt * node];

                        /* Break if time is past next unit of time */
                        if (isinf(sa.t_time[node]) || sa.t_time[node] >= sa.next_unit_of_time) {
                            sa.t_time[node] = sa.next_unit_of_time;
                            break;
                        }

                        /* 1b) Determine which transitions that occur */
                        tr = ma.reactNode[sa.Nt * node]%sa.Nt;

                        /* 1c) Update the state of the node */
                        for (j = sa.jcS[tr]; j < sa.jcS[tr + 1]; j++) {
                            sa.u[node * sa.Nc + sa.irS[j]] += sa.prS[j];
                            if (sa.u[node * sa.Nc + sa.irS[j]] < 0) {
                                SimInf_print_status(sa.Nc, &sa.u[node * sa.Nc],
                                                    sa.Nd, &sa.v[node * sa.Nd],
                                                    sa.Nld, &sa.ldata[node * sa.Nld],
                                                    sa.Ni + node, sa.t_time[node],
                                                    0, tr);
                                sa.error = SIMINF_ERR_NEGATIVE_STATE;
                            }
                        }


                        /* 1d) update dependent transitions events. */
                        for (ii = sa.jcG[tr]; ii < sa.jcG[tr + 1]; ii++){
                            j = sa.irG[ii];
                            if (j != tr) { /*see code underneath */
                                old_t_rate = sa.t_rate[node * sa.Nt + j];
                                /* const double rate */
                                rate = (*sa.tr_fun[j])(
                                    &sa.u[node * sa.Nc], &sa.v[node * sa.Nd],
                                    &sa.ldata[node * sa.Nld], sa.gdata,
                                    sa.t_time[node]);

                                sa.t_rate[node * sa.Nt + j] = rate;

                                if (!R_FINITE(rate) || rate < 0.0) {
                                    SimInf_print_status(sa.Nc, &sa.u[node * sa.Nc],
                                                        sa.Nd, &sa.v[node * sa.Nd],
                                                        sa.Nld, &sa.ldata[node * sa.Nld],
                                                        sa.Ni + node, sa.t_time[node],
                                                        rate, j);
                                    sa.error = SIMINF_ERR_INVALID_RATE;
                                }

                                /* update times and reorder the heap */
                                calcTimes(&ma.reactTimes[sa.Nt * node + ma.reactHeap[sa.Nt * node + j]],
                                          &ma.reactInf[sa.Nt * node + j],
                                          sa.t_time[node],
                                          old_t_rate,
                                          sa.t_rate[node * sa.Nt + j],
                                          ma.rng_vec[sa.Nt * node + j]);
                                update(ma.reactHeap[sa.Nt * node + j], &ma.reactTimes[sa.Nt * node],
                                       &ma.reactNode[sa.Nt * node], &ma.reactHeap[sa.Nt * node], ma.reactHeapSize);
                            }
                        }
                        /* finish with j = re (the one that just happened), which need
                           not be in the dependency graph but must be updated  nevertheless */
                        j = tr;
                        old_t_rate = sa.t_rate[node * sa.Nt + j];
                        rate = (*sa.tr_fun[j])(&sa.u[node * sa.Nc], &sa.v[node * sa.Nd],
                                               &sa.ldata[node * sa.Nld], sa.gdata,
                                               sa.t_time[node]);
                        sa.t_rate[node * sa.Nt + j] = rate;

                        if (!R_FINITE(rate) || rate < 0.0) {
                            SimInf_print_status(sa.Nc, &sa.u[node * sa.Nc],
                                                sa.Nd, &sa.v[node * sa.Nd],
                                                sa.Nld, &sa.ldata[node * sa.Nld],
                                                sa.Ni + node, sa.t_time[node],
                                                rate, j);
                            sa.error = SIMINF_ERR_INVALID_RATE;
                        }

                        /* update times and reorder the heap */
                        calcTimes(&ma.reactTimes[sa.Nt * node + ma.reactHeap[sa.Nt * node + j]],
                                  &ma.reactInf[sa.Nt * node + j],
                                  sa.t_time[node],
                                  old_t_rate,
                                  sa.t_rate[node * sa.Nt + j],
                                  ma.rng_vec[sa.Nt * node + j]);
                        update(ma.reactHeap[sa.Nt * node + j], &ma.reactTimes[sa.Nt * node],
                               &ma.reactNode[sa.Nt * node], &ma.reactHeap[sa.Nt * node], ma.reactHeapSize);

                    }
                }

                *&model[i] = sa;
                *&method[i] = ma;

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
                SimInf_compartment_model sa = *&model[i];
                SimInf_aem_arguments ma = *&method[i];

                /* (4) Incorporate model specific actions after each
                 * timestep e.g. update the infectious pressure
                 * variable. Moreover, update transition rates in
                 * nodes that are indicated for update */
                for (node = 0; node < sa.Nn; node++) {
                    const int rc = sa.pts_fun(
                        &sa.v_new[node * sa.Nd], &sa.u[node * sa.Nc],
                        &sa.v[node * sa.Nd], &sa.ldata[node * sa.Nld],
                        sa.gdata, sa.Ni + node, sa.tt);

                    if (rc < 0) {
                        sa.error = rc;
                        break;
                    } else if (rc > 0 || sa.update_node[node]) {
                        /* Update transition rates */
                        int j = 0;
                        for (; j < sa.Nt; j++) {
                            const double old = sa.t_rate[node * sa.Nt + j];
                            const double rate = (*sa.tr_fun[j])(
                                &sa.u[node * sa.Nc], &sa.v_new[node * sa.Nd],
                                &sa.ldata[node * sa.Nld], sa.gdata, sa.tt);

                            sa.t_rate[node * sa.Nt + j] = rate;

                            if (!R_FINITE(rate) || rate < 0.0) {
                                SimInf_print_status(sa.Nc, &sa.u[node * sa.Nc],
                                                    sa.Nd, &sa.v[node * sa.Nd],
                                                    sa.Nld, &sa.ldata[node * sa.Nld],
                                                    sa.Ni + node, sa.tt, rate, j);
                                sa.error = SIMINF_ERR_INVALID_RATE;
                            }

			    /* Update times and reorder heap */
			    calcTimes(&ma.reactTimes[sa.Nt * node + ma.reactHeap[sa.Nt * node + j]],
				      &ma.reactInf[sa.Nt * node + j],
				      sa.t_time[node],
				      old,
				      sa.t_rate[node * sa.Nt + j],
				      ma.rng_vec[sa.Nt * node + j]);

			    update(ma.reactHeap[sa.Nt * node + j], &ma.reactTimes[sa.Nt * node],
                                   &ma.reactNode[sa.Nt * node], &ma.reactHeap[sa.Nt * node], ma.reactHeapSize);
                        }

                        sa.update_node[node] = 0;
                    }
                }

                /* (5) The global time now equals next unit of time. */
                sa.tt = sa.next_unit_of_time;
                sa.next_unit_of_time += 1.0;

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
                while (sa.U && sa.U_it < sa.tlen && sa.tt > sa.tspan[sa.U_it])
                    memcpy(&sa.U[sa.Nc * ((sa.Ntot * sa.U_it++) + sa.Ni)],
                           sa.u, sa.Nn * sa.Nc * sizeof(int));
                /* Copy continuous state to V */
                while (sa.V && sa.V_it < sa.tlen && sa.tt > sa.tspan[sa.V_it])
                    memcpy(&sa.V[sa.Nd * ((sa.Ntot * sa.V_it++) + sa.Ni)],
                           sa.v_new, sa.Nn * sa.Nd * sizeof(double));

                *&model[i] = sa;
                *&method[i] = ma;
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
 * Free allocated memory for an epidemiological compartment
 * model.
 *
 * @param method the data structure to free
 * @param model structure with data about the model
 * @param Nthread number of threads that was used during simulation.
 */
static void SimInf_aem_arguments_free(
    SimInf_aem_arguments *method, SimInf_compartment_model *model ,int Nthread)
{
    if (method) {
        int i;

        for (i = 0; i < Nthread; i++) {
            SimInf_aem_arguments *m = &method[i];
            SimInf_compartment_model *mod = &model[i];

            if (m) {
                /* AEM variables */
                if (m->rng_vec) {
                    int j;
                    for (j = 0; j < mod->Nn * mod->Nt; j++)
                        gsl_rng_free(m->rng_vec[j]);
                }
                m->rng_vec = NULL;
                free(m->reactHeap);
                m->reactHeap = NULL;
                free(m->reactInf);
                m->reactInf = NULL;
                free(m->reactNode);
                m->reactNode = NULL;
                free(m->reactTimes);
                m->reactTimes = NULL;
            }
        }
        free(method);
    }
}

/**
 * Create and initialize data for an epidemiological compartment
 * model. The generated model must be freed by the user.
 *
 * @param out the resulting data structure.
 * @param model structure with data about the model
 * @param Nthread the number of threads available
 * @param rng random number generator.
 * @return 0 or SIMINF_ERR_ALLOC_MEMORY_BUFFER
 */
static int SimInf_aem_arguments_create(
    SimInf_aem_arguments **out, SimInf_compartment_model *model, int Nthread, gsl_rng *rng)
{
    int i;
    SimInf_aem_arguments *method = NULL;

    method = calloc(Nthread, sizeof(SimInf_aem_arguments));
    if(!method)
        goto on_error; /* #nocov */

    for (i = 0; i < Nthread; i++) {
        int node;
        SimInf_compartment_model *m = &model[i];
        /* Binary heap storing all reaction events */
        /* we have one for each node. Heap is thus only the size of the # transitions */
        method[i].reactHeapSize = m->Nt;
        method[i].reactNode = malloc(m->Nn * m->Nt * sizeof(int));
        if (!method[i].reactNode)
            goto on_error; /* #nocov */

        method[i].reactHeap = malloc(m->Nn * m->Nt * sizeof(int));
        if (!method[i].reactHeap)
            goto on_error; /* #nocov */

        method[i].reactTimes = malloc(m->Nn * m->Nt * sizeof(double));
        if (!method[i].reactTimes)
            goto on_error; /* #nocov */

        method[i].reactInf = calloc(m->Nn * m->Nt, sizeof(double));
        if (!method[i].reactInf)
            goto on_error; /* #nocov */

        /* random generator for sample select with 1 per transition in each node */
        method[i].rng_vec = malloc(m->Nn * m->Nt * sizeof(gsl_rng*));
        if (!method[i].rng_vec)
            goto on_error; /* #nocov */

        for (node = 0; node < m->Nn; node++) {
            int trans;
            for (trans = 0; trans < m->Nt; trans++) {
                /* Random number generator */
                method[i].rng_vec[m->Nt * node + trans] = gsl_rng_alloc(gsl_rng_mt19937);
                if (!method[i].rng_vec[m->Nt * node + trans])
                    goto on_error; /* #nocov */

                gsl_rng_set(method[i].rng_vec[m->Nt * node + trans],
                            gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
            }
        }
    }

    *out = method;

    return 0;

on_error:                                              /* #nocov */
    SimInf_aem_arguments_free(method, model, Nthread); /* #nocov */
    return SIMINF_ERR_ALLOC_MEMORY_BUFFER;             /* #nocov */
}

/**
 * Initialize and run siminf solver
 *
 * @param args Structure with data for the solver.
 * @return 0 if Ok, else error code.
 */
int attribute_hidden SimInf_run_solver_aem(SimInf_solver_args *args)
{
    int error = 0;
    gsl_rng *rng = NULL;
    SimInf_scheduled_events *events = NULL;
    SimInf_compartment_model *model = NULL;
    SimInf_aem_arguments *method = NULL;

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

    error = SimInf_aem_arguments_create(&method, model, args->Nthread, rng);
    if (error)
        goto cleanup; /* #nocov */

    error = SimInf_solver_aem(model, method, events, args->Nthread);

cleanup:
    gsl_rng_free(rng);
    SimInf_scheduled_events_free(events);
    SimInf_aem_arguments_free(method, model, args->Nthread);
    SimInf_compartment_model_free(model);

    return error;
}
