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
#include "SimInf_solver_aem.h"
#include "binheap.h"

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
void calcTimes(double* time, double* infTime, double tt, double old_rate,
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
    SimInf_compartment_model *sim_args, SimInf_model_events *events,
    int *uu, int *update_node, int Nthread)
{
    int k;

    #pragma omp parallel
    {
        int i;

        #pragma omp for
        for (i = 0; i < Nthread; i++) {
            int node;
            SimInf_compartment_model sa = *&sim_args[i];

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

                    if (!isfinite(rate) || rate < 0.0)
                        sa.errcode = SIMINF_ERR_INVALID_RATE;

                    /* calculate time until next transition j event */
                    sa.reactTimes[sa.Nt*node+j] =  -log(gsl_rng_uniform_pos(sa.rng_vec[sa.Nt*node+j]))/rate + sa.tt;
                    if (sa.reactTimes[sa.Nt*node+j] <= 0.0)
                        sa.reactTimes[sa.Nt*node+j] = INFINITY;

                    sa.reactHeap[sa.Nt*node+j] = sa.reactNode[sa.Nt*node+j] = j;
                }

                /* Initialize reaction heap */
                initialize_heap(&sa.reactTimes[sa.Nt*node], &sa.reactNode[sa.Nt*node],
                                &sa.reactHeap[sa.Nt*node], sa.reactHeapSize);
                sa.t_time[node] = sa.tt;
	    }

	    *&sim_args[i] = sa;
        }
    }

    /* Check for error during initialization. */
    for (k = 0; k < Nthread; k++)
        if (sim_args[k].errcode)
            return sim_args[k].errcode;

    /* Main loop. */
    for (;;) {
        #pragma omp parallel
        {
            int i;

            #pragma omp for
            for (i = 0; i < Nthread; i++) {
                int node;
                SimInf_compartment_model sa = *&sim_args[i];

                /* (1) Handle internal epidemiological model,
                 * continuous-time Markov chain. */
                for (node = 0; node < sa.Nn && !sa.errcode; node++) {
                    for (;;) {
                        int ii,j,tr;
                        double old_t_rate,rate;

                        /* 1a) Step time forward until next event */
                        sa.t_time[node] = sa.reactTimes[sa.Nt * node];

                        /* If time is past next day break */
                        if (isinf(sa.t_time[node]) || sa.t_time[node] >= sa.next_day) {
                            sa.t_time[node] = sa.next_day;
                            break;
                        }

                        /* 1b) Determine which transitions that occur */
                        tr = sa.reactNode[sa.Nt * node]%sa.Nt;

                        /* 1c) Update the state of the node */
                        for (j = sa.jcS[tr]; j < sa.jcS[tr + 1]; j++) {
                            sa.u[node * sa.Nc + sa.irS[j]] += sa.prS[j];
                            if (sa.u[node * sa.Nc + sa.irS[j]] < 0)
                                sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
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

                                if (!isfinite(rate) || rate < 0.0)
                                    sa.errcode = SIMINF_ERR_INVALID_RATE;
                                /* update times and reorder the heap */
                                calcTimes(&sa.reactTimes[sa.Nt * node + sa.reactHeap[sa.Nt * node + j]],
                                          &sa.reactInf[sa.Nt * node + j],
                                          sa.t_time[node],
                                          old_t_rate,
                                          sa.t_rate[node * sa.Nt + j],
                                          sa.rng_vec[sa.Nt * node + j]);
                                update(sa.reactHeap[sa.Nt * node + j], &sa.reactTimes[sa.Nt * node],
                                       &sa.reactNode[sa.Nt * node], &sa.reactHeap[sa.Nt * node], sa.reactHeapSize);
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

                        if (!isfinite(rate) || rate < 0.0)
                            sa.errcode = SIMINF_ERR_INVALID_RATE;

                        /* update times and reorder the heap */
                        calcTimes(&sa.reactTimes[sa.Nt * node + sa.reactHeap[sa.Nt * node + j]],
                                  &sa.reactInf[sa.Nt * node + j],
                                  sa.t_time[node],
                                  old_t_rate,
                                  sa.t_rate[node * sa.Nt + j],
                                  sa.rng_vec[sa.Nt * node + j]);
                        update(sa.reactHeap[sa.Nt * node + j], &sa.reactTimes[sa.Nt * node],
                               &sa.reactNode[sa.Nt * node], &sa.reactHeap[sa.Nt * node], sa.reactHeapSize);

                    }
                }

                *&sim_args[i] = sa;

                /* (2) Incorporate all scheduled E1 events */
                SimInf_process_E1_events(&sim_args[i], &events[i], uu, update_node);
	    }

            #pragma omp barrier

            #pragma omp master
            {
                /* (3) Incorporate all scheduled E2 events */
                SimInf_process_E2_events(sim_args, events, uu, update_node);
            }

            #pragma omp barrier

            #pragma omp for
            for (i = 0; i < Nthread; i++) {
                int node;
                SimInf_compartment_model sa = *&sim_args[i];

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
                        sa.errcode = rc;
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

                            if (!isfinite(rate) || rate < 0.0)
                                sa.errcode = SIMINF_ERR_INVALID_RATE;

			    /* Update times and reorder heap */
			    calcTimes(&sa.reactTimes[sa.Nt * node + sa.reactHeap[sa.Nt * node + j]],
				      &sa.reactInf[sa.Nt * node + j],
				      sa.t_time[node],
				      old,
				      sa.t_rate[node * sa.Nt + j],
				      sa.rng_vec[sa.Nt * node + j]);

			    update(sa.reactHeap[sa.Nt * node + j], &sa.reactTimes[sa.Nt * node],
                                   &sa.reactNode[sa.Nt * node], &sa.reactHeap[sa.Nt * node], sa.reactHeapSize);
                        }

                        sa.update_node[node] = 0;
                    }
                }

                /* (5) The global time now equals next_day. */
                sa.tt = sa.next_day;
                sa.next_day += 1.0;

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

                *&sim_args[i] = sa;
            }
        }

        /* 6b) Handle the case where the solution is stored in a sparse
         * matrix */
        SimInf_store_solution_sparse(sim_args);

        /* Swap the pointers to the continuous state variable so that
         * 'v' equals 'v_new'. Moreover, check for error. */
        for (k = 0; k < Nthread; k++) {
            double *v_tmp = sim_args[k].v;
            sim_args[k].v = sim_args[k].v_new;
            sim_args[k].v_new = v_tmp;
            if (sim_args[k].errcode)
                return sim_args[k].errcode;
        }

        /* If the simulation has reached the final time, exit. */
        if (sim_args[0].U_it >= sim_args[0].tlen)
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
int SimInf_run_solver_aem(SimInf_solver_args *args)
{
    int i, error;
    gsl_rng *rng = NULL;
    SimInf_model_events *events = NULL;
    SimInf_compartment_model *sim_args = NULL;
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
        &sim_args, args, uu, vv_1, vv_2, update_node);
    if (error)
        goto cleanup;

    error = SimInf_model_events_create(&events, args, rng);
    if (error)
        goto cleanup;

    for (i = 0; i < args->Nthread; i++) {
        int node;

        /* start AEM specific */
        /* Binary heap storing all reaction events */
        /* we have one for each node. Heap is thus only the size of the # transitions */
        sim_args[i].reactHeapSize = args->Nt;
        sim_args[i].reactNode = malloc(args->Nn * args->Nt * sizeof(int));
        if (!sim_args[i].reactNode) {
            error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        sim_args[i].reactHeap = malloc(args->Nn * args->Nt * sizeof(int));
        if (!sim_args[i].reactHeap) {
            error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        sim_args[i].reactTimes = malloc(args->Nn * args->Nt * sizeof(double));
        if (!sim_args[i].reactTimes) {
            error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        sim_args[i].reactInf = calloc(args->Nn * args->Nt, sizeof(double));
        if (!sim_args[i].reactInf) {
            error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        /* random generator for sample select with 1 per transition in each node */
        sim_args[i].rng_vec = malloc(args->Nn * args->Nt * sizeof(gsl_rng*));
        if (!sim_args[i].rng_vec) {
            error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        for (node = 0; node < args->Nn; node++) {
            int trans;
            for (trans = 0; trans < sim_args->Nt; trans++) {
                /* Random number generator */
                sim_args[i].rng_vec[args->Nt * node + trans] = gsl_rng_alloc(gsl_rng_mt19937);
                if (!sim_args[i].rng_vec[args->Nt * node + trans]) {
                    error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
                    goto cleanup;
                }
                if (!sim_args[i].rng_vec[args->Nt * node + trans]) {
                    error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
                    goto cleanup;
                }
                gsl_rng_set(sim_args[i].rng_vec[args->Nt * node + trans],
                            gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
            }
        }
        /* end AEM specific */

    }
    if (error)
        goto cleanup;

    error = SimInf_solver_aem(sim_args, events, uu, update_node, args->Nthread);

cleanup:
    if (uu) {
        free(uu);
        uu = NULL;
    }

    if (vv_1) {
        free(vv_1);
        vv_1 = NULL;
    }

    if (vv_2) {
        free(vv_2);
        vv_2 = NULL;
    }

    if (update_node) {
        free(update_node);
        update_node = NULL;
    }

    if (rng)
        gsl_rng_free(rng);

    if (events) {
        for (i = 0; i < args->Nthread; i++)
            SimInf_free_model_events(&events[i]);
        free(events);
        events = NULL;
    }

    if (sim_args) {
        for (i = 0; i < args->Nthread; i++)
            SimInf_free_args(&sim_args[i]);
        free(sim_args);
        sim_args = NULL;
    }

    return error;
}
