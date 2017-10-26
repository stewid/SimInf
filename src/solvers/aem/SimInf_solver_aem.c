/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2017  Robin Eriksson
 *  Copyright (C) 2015 - 2017 Stefan Engblom
 *  Copyright (C) 2015 - 2017 Stefan Widgren
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

        if (infTime[0] == 0) { // Waking up first time
            time[0] = -log(1.0 - gsl_rng_uniform_pos(rng)) / new_rate + tt;
        } else { // Waking up the 2nd..nth time
            if (new_rate > 0.0) {
                time[0] = tt + (infTime[0] / new_rate);
            }
        }
    }
    else {
        if (new_rate >= DBL_MIN) {
            if (time[0] == tt) // Regular update of current event
                time[0] = -log(1.0 - gsl_rng_uniform_pos(rng)) / new_rate + tt;
            else
                // Regular update of dependent events (rescaling)
                time[0] = ((old_rate / new_rate) * (time[0] - tt)) + tt;
        } else { // Next event time set to infinity

            infTime[0] = (oldtime - tt) * old_rate;
            time[0] = INFINITY;
        }
    }
}

/**
 * Siminf solver
 *
 * @return 0 if Ok, else error code.
 */
static int SimInf_solver_aem(
        SimInf_thread_args *sim_args, int *uu, int *update_node, int Nthread)
{
    int k;

    #pragma omp parallel
    {
        int i;

        #pragma omp for
        for (i = 0; i < Nthread; i++) {
            int node;
            SimInf_thread_args sa = *&sim_args[i];

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

                    if(!isfinite(rate) || rate < 0.0)
                        sa.errcode = SIMINF_ERR_INVALID_RATE;

                    /* calculate time until next transition j event */
                    sa.reactTimes[sa.Nt*node+j] =  -log(1.0-gsl_rng_uniform_pos(sa.rng_vec[sa.Nt*node+j]))/rate + sa.tt;
                    if(sa.reactTimes[sa.Nt*node+j] <= 0.0)
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
                SimInf_thread_args sa = *&sim_args[i];
                SimInf_scheduled_events e1 = *sa.E1;

                /* (1) Handle internal epidemiological model,
                 * continuous-time Markov chain. */
                for (node = 0; node < sa.Nn && !sa.errcode; node++) {
                    for (;;) {
                        int i,j,tr; /*tr <- re in urdme:aem.c */
                        double old_t_rate,rate;



                        /* AEM solver */

                        /* 1a) Step time forward until next event */
                        sa.t_time[node] = sa.reactTimes[sa.Nt*node];

                        /* If time is past next day break */
                        if(sa.t_time[node] >= sa.next_day || isinf(sa.t_time[node])){
                            sa.t_time[node] = sa.next_day;
                            break;
                        }

                        /* 1b) Determine which transistion that occur */
                        tr = sa.reactNode[sa.Nt*node]%sa.Nt;

                        /* 1c) Update the state of the node */
                        for (j = sa.jcS[tr]; j < sa.jcS[tr + 1]; j++) {
                            sa.u[node * sa.Nc + sa.irS[j]] += sa.prS[j];
                            if (sa.u[node * sa.Nc + sa.irS[j]] < 0)
                                sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                        }

                        /* 1d) update dependent transitions events. */
                        for (i = sa.jcG[tr]; i < sa.jcG[tr + 1]; i++) {
                            j = sa.irG[i];
                            if(j != tr) { /*see code underneath */
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
                                calcTimes(&sa.reactTimes[sa.Nt*node + sa.reactHeap[sa.Nt*node+j]],
                                          &sa.reactInf[sa.Nt*node+j],
                                          sa.t_time[node],
                                          old_t_rate,
                                          sa.t_rate[node * sa.Nt + j],
                                          sa.rng_vec[sa.Nt*node+j]);
                                update(sa.reactHeap[sa.Nt*node+j], &sa.reactTimes[sa.Nt*node],
                                       &sa.reactNode[sa.Nt*node], &sa.reactHeap[sa.Nt*node], sa.reactHeapSize);
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
                        calcTimes(&sa.reactTimes[sa.Nt*node + sa.reactHeap[sa.Nt*node+j]],
                                  &sa.reactInf[sa.Nt*node+j],
                                  sa.t_time[node],
                                  old_t_rate,
                                  sa.t_rate[node * sa.Nt + j],
                                  sa.rng_vec[sa.Nt*node+j]);
                        update(sa.reactHeap[sa.Nt*node+j], &sa.reactTimes[sa.Nt*node],
                               &sa.reactNode[sa.Nt*node], &sa.reactHeap[sa.Nt*node], sa.reactHeapSize);
                    }
		}

                /* (2) Incorporate all scheduled E1 events */
                while (sa.E1_index < e1.len &&
                       sa.tt >= e1.time[sa.E1_index] &&
                       !sa.errcode)
                {
                    const int j = sa.E1_index;
                    const int s = e1.select[j];

                    if (e1.event[j] == ENTER_EVENT) {
                        /* All individuals enter first non-zero
                         * compartment, i.e. a non-zero entry in
                         * element in the select column. */
                        if (sa.jcE[s] < sa.jcE[s + 1]) {
                            uu[e1.node[j] * sa.Nc + sa.irE[sa.jcE[s]]] += e1.n[j];
                            if (uu[e1.node[j] * sa.Nc + sa.irE[sa.jcE[s]]] < 0)
                                sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                        }
                    } else {
                        sa.errcode = SimInf_sample_select(
                            sa.irE, sa.jcE, sa.Nc, uu, e1.node[j],
                            e1.select[j], e1.n[j], e1.proportion[j],
                            sa.individuals, sa.u_tmp, sa.rng);

                        if (sa.errcode)
                            break;

                        if (e1.event[j] == EXIT_EVENT) {
                            int ii;

                            for (ii = sa.jcE[s]; ii < sa.jcE[s + 1]; ii++) {
                                const int jj = sa.irE[ii];
                                const int kk = e1.node[j] * sa.Nc + jj;

                                /* Remove individuals from node */
                                uu[kk] -= sa.individuals[jj];
                                if (uu[kk] < 0) {
                                    sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                                    break;
                                }
                            }
                        } else { /* INTERNAL_TRANSFER_EVENT */
                            int ii;

                            for (ii = sa.jcE[s]; ii < sa.jcE[s + 1]; ii++) {
                                const int jj = sa.irE[ii];
                                const int kk = e1.node[j] * sa.Nc + jj;
                                const int ll = sa.N[e1.shift[j] * sa.Nc + jj];

                                /* Add individuals to new compartments
                                 * in node */
                                uu[kk + ll] += sa.individuals[jj];
                                if (uu[kk + ll] < 0) {
                                    sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                                    break;
                                }

                                /* Remove individuals from previous
                                 * compartments in node */
                                uu[kk] -= sa.individuals[jj];
                                if (uu[kk] < 0) {
                                    sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                                    break;
                                }
                            }
                        }
                    }

                    /* Indicate node for update */
                    update_node[e1.node[j]] = 1;
                    sa.E1_index++;
                }

                *&sim_args[i] = sa;
	    }


            #pragma omp barrier

            #pragma omp master
            {
                SimInf_thread_args sa = *&sim_args[0];
                SimInf_scheduled_events e2 = *sa.E2;

                /* (3) Incorporate all scheduled E2 events */
                while (sa.E2_index < e2.len &&
                       sa.tt >= e2.time[sa.E2_index] &&
                       !sa.errcode)
                {
                    sa.errcode = SimInf_sample_select(
                        sa.irE, sa.jcE, sa.Nc, uu, e2.node[sa.E2_index],
                        e2.select[sa.E2_index], e2.n[sa.E2_index],
                        e2.proportion[sa.E2_index], sa.individuals,
                        sa.u_tmp, sa.rng);

                    if (sa.errcode)
                        break;

                    for (i = sa.jcE[e2.select[sa.E2_index]];
                         i < sa.jcE[e2.select[sa.E2_index] + 1];
                         i++)
                    {
                        const int jj = sa.irE[i];
                        const int k1 = e2.dest[sa.E2_index] * sa.Nc + jj;
                        const int k2 = e2.node[sa.E2_index] * sa.Nc + jj;
                        const int ll = e2.shift[sa.E2_index] < 0 ? 0 :
                            sa.N[e2.shift[sa.E2_index] * sa.Nc + jj];

                        /* Add individuals to dest */
                        uu[k1 + ll] += sa.individuals[jj];
                        if (uu[k1 + ll] < 0) {
                            sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                            break;
                        }

                        /* Remove individuals from node */
                        uu[k2] -= sa.individuals[jj];
                        if (uu[k2] < 0) {
                            sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                            break;
                        }
                    }

                    /* Indicate node and dest for update */
                    update_node[e2.node[sa.E2_index]] = 1;
                    update_node[e2.dest[sa.E2_index]] = 1;
                    sa.E2_index++;
                }

                *&sim_args[0] = sa;
            }

            #pragma omp barrier

            #pragma omp for
            for (i = 0; i < Nthread; i++) {
                int node;
                SimInf_thread_args sa = *&sim_args[i];

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
			    calcTimes(&sa.reactTimes[sa.Nt*node + sa.reactHeap[sa.Nt*node+j]],
				      &sa.reactInf[sa.Nt*node+j],
				      sa.t_time[node],
				      old,
				      sa.t_rate[node * sa.Nt + j],
				      sa.rng_vec[sa.Nt*node+j]);

			    update(sa.reactHeap[sa.Nt*node+j], &sa.reactTimes[sa.Nt*node],
                                   &sa.reactNode[sa.Nt*node], &sa.reactHeap[sa.Nt*node], sa.reactHeapSize);
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
        while (!sim_args[0].U && sim_args[0].U_it < sim_args[0].tlen &&
               sim_args[0].tt > sim_args[0].tspan[sim_args[0].U_it]) {
            int j;

            /* Copy compartment state to U_sparse */
            for (j = sim_args[0].jcU[sim_args[0].U_it];
                 j < sim_args[0].jcU[sim_args[0].U_it + 1]; j++)
                sim_args[0].prU[j] = sim_args[0].u[sim_args[0].irU[j]];
            sim_args[0].U_it++;
        }

        while (!sim_args[0].V && sim_args[0].V_it < sim_args[0].tlen &&
               sim_args[0].tt > sim_args[0].tspan[sim_args[0].V_it]) {
            int j;

            /* Copy continuous state to V_sparse */
            for (j = sim_args[0].jcV[sim_args[0].V_it];
                 j < sim_args[0].jcV[sim_args[0].V_it + 1]; j++)
                sim_args[0].prV[j] = sim_args[0].v_new[sim_args[0].irV[j]];
            sim_args[0].V_it++;
        }

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
    int i, errcode;
    gsl_rng *rng = NULL;
    SimInf_thread_args *sim_args = NULL;
    int *uu = NULL, *update_node = NULL;
    double *vv_1 = NULL, *vv_2 = NULL;

    /* Set compartment state to the initial state. */
    uu = malloc(args->Nn * args->Nc * sizeof(int));
    if (!uu) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
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
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
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
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!rng) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    gsl_rng_set(rng, args->seed);

    sim_args = calloc(args->Nthread, sizeof(SimInf_thread_args));
    if (!sim_args) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }


    for (i = 0; i < args->Nthread; i++) {
        /* Random number generator */
        sim_args[i].rng = gsl_rng_alloc(gsl_rng_mt19937);
        if (!sim_args[i].rng) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }
        gsl_rng_set(sim_args[i].rng, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));

        /* Constants */
        sim_args[i].Ntot = args->Nn;
        sim_args[i].Ni = i * (args->Nn / args->Nthread);
        sim_args[i].Nn = args->Nn / args->Nthread;
        if (i == (args->Nthread - 1))
            sim_args[i].Nn += (args->Nn % args->Nthread);
        sim_args[i].Nt = args->Nt;
        sim_args[i].Nc = args->Nc;
        sim_args[i].Nd = args->Nd;
        sim_args[i].Nld = args->Nld;

        /* Sparse matrices */
        sim_args[i].irG = args->irG;
        sim_args[i].jcG = args->jcG;
        sim_args[i].irS = args->irS;
        sim_args[i].jcS = args->jcS;
        sim_args[i].prS = args->prS;
        sim_args[i].irE = args->irE;
        sim_args[i].jcE = args->jcE;

        /* Callbacks */
        sim_args[i].tr_fun = args->tr_fun;
        sim_args[i].pts_fun = args->pts_fun;

        /* Keep track of time */
        sim_args[i].tt = args->tspan[0];
        sim_args[i].next_day = floor(sim_args[i].tt) + 1.0;
        sim_args[i].tspan = args->tspan;
        sim_args[i].tlen = args->tlen;
        sim_args[i].U_it = 1;
        sim_args[i].V_it = 1;

        /* start AEM specific
         */
        /* Binary heap storing all reaction events */
        /* we have one for each node. Heap is thus only the size of the # transitions */
        sim_args[i].reactHeapSize = sim_args[i].Nt;

        sim_args[i].reactNode = (int *)malloc(sim_args[i].Nn*sim_args[i].Nt*sizeof(int));
        sim_args[i].reactHeap = (int *)malloc(sim_args[i].Nn*sim_args[i].Nt*sizeof(int));
        sim_args[i].reactTimes = (double *)malloc(sim_args[i].Nn*sim_args[i].Nt*sizeof(double));
        /* Initialize with zero elements */
        sim_args[i].reactInf = (double *)calloc(sim_args[i].Nn*sim_args[i].Nt,sizeof(double));

        /* random generator for sample select with 1 per transition in each node */
        sim_args[i].rng_vec = (gsl_rng **)malloc(sim_args[i].Nn*sim_args[i].Nt*sizeof(gsl_rng*));

        for(int node = 0; node < sim_args[i].Nn; node++){
            /* sim_args[i].rng_vec[node] = (gsl_rng **)malloc(sim_args[i].Nt*sizeof(gsl_rng*)); */
            for(int trans = 0; trans < sim_args[i].Nt; trans++){
                /* Random number generator */
                sim_args[i].rng_vec[sim_args[i].Nt*node + trans] = gsl_rng_alloc(gsl_rng_mt19937);
                if (!sim_args[i].rng_vec[sim_args[i].Nt*node + trans]) {
                    errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
                    goto cleanup;
                }
                gsl_rng_set(sim_args[i].rng_vec[sim_args[i].Nt*node + trans],
                            gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
            }
        }
        /* end AEM specific
         */

        /* Data vectors */
        sim_args[i].N = args->N;
        if (args->U) {
            sim_args[i].U = args->U;
        } else if (i == 0) {
            sim_args[i].irU = args->irU;
            sim_args[i].jcU = args->jcU;
            sim_args[i].prU = args->prU;
        }
        sim_args[i].u = &uu[sim_args[i].Ni * args->Nc];
        if (args->V) {
            sim_args[i].V = args->V;
        } else if (i == 0) {
            sim_args[i].irV = args->irV;
            sim_args[i].jcV = args->jcV;
            sim_args[i].prV = args->prV;
        }
        sim_args[i].v = &vv_1[sim_args[i].Ni * args->Nd];
        sim_args[i].v_new = &vv_2[sim_args[i].Ni * args->Nd];
        sim_args[i].ldata = &(args->ldata[sim_args[i].Ni * sim_args[i].Nld]);
        sim_args[i].gdata = args->gdata;
        sim_args[i].update_node = &update_node[sim_args[i].Ni];

        /* Scheduled events */
        sim_args[i].E1 = calloc(1, sizeof(SimInf_scheduled_events));
        if (!sim_args[i].E1) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        if (i == 0) {
            sim_args[i].E2 = calloc(1, sizeof(SimInf_scheduled_events));
            if (!sim_args[i].E2) {
                errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
                goto cleanup;
            }
        }

        sim_args[i].individuals = calloc(args->Nc, sizeof(int));
        if (!sim_args[i].individuals) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        sim_args[i].u_tmp = calloc(args->Nc, sizeof(int));
        if (!sim_args[i].u_tmp) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

	/* Create transition rate matrix (Nt X Nn) and total rate
         * vector. In t_rate we store all propensities for state
         * transitions, and in sum_t_rate the sum of propensities
         * in every node. */
        sim_args[i].t_rate = malloc(args->Nt * sim_args[i].Nn * sizeof(double));
        if (!sim_args[i].t_rate) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }
        sim_args[i].sum_t_rate = malloc(sim_args[i].Nn * sizeof(double));
        if (!sim_args[i].sum_t_rate) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }
        sim_args[i].t_time = malloc(sim_args[i].Nn * sizeof(double));
        if (!sim_args[i].t_time) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }
    }

    /* Split scheduled events into E1 and E2 events. */
    errcode = SimInf_split_events(sim_args,
        args->len, args->event, args->time, args->node, args->dest, args->n,
        args->proportion, args->select, args->shift, args->Nn, args->Nthread);
    if (errcode)
        goto cleanup;

    errcode = SimInf_solver_aem(sim_args, uu, update_node, args->Nthread);

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

    if (sim_args) {
        for (i = 0; i < args->Nthread; i++)
            SimInf_free_args(&sim_args[i]);
        free(sim_args);
        sim_args = NULL;
    }

    return errcode;
}
