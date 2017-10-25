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
 * Event types
 *
 * EXIT_EVENT (0): Exit events are events that remove individuals from
 * a node.
 *
 * ENTER_EVENT (1): Enter events are events that introduce new
 * individuals into a node. All individuals enter first non-zero
 * compartment, i.e. a non-zero entry in element in the select column.
 *
 * INTERNAL_TRANSFER_EVENT (2): Internal transfer events are events
 * that change the number of individuals in the compartments whithin
 * one node e.g. aging of n individuals from age_1 to age_2 in a model
 * with age categories.
 *
 * EXTERNAL_TRANSFER_EVENT (3): External transfer events are events
 * that move individuals from compartments in one node to compartments
 * in another node e.g. moving n individuals from node A to node B.
 */
enum {EXIT_EVENT,
      ENTER_EVENT,
      INTERNAL_TRANSFER_EVENT,
      EXTERNAL_TRANSFER_EVENT};

/**
 * Structure to hold thread specific data/arguments for simulation.
 */
typedef struct SimInf_thread_args
{
    /*** Constants ***/
    int Ntot;  /**< Total number of nodes. */
    int Ni;    /**< Index to first node in thread in the global set of
                *  of nodes. */
    int Nn;    /**< Number of nodes in thread. */
    int Nt;    /**< Total number of different transitions. */
    int Nc;    /**< Number of compartments in each node. */
    int Nd;    /**< Number of continuous state variables. */
    int Nld;   /**< Length of the local data vector 'ldata' for each
                *   node. The 'ldata' vector is sent to propensities
                *   and the post time step function. */

    /*** Sparse matrices ***/
    const int *irG; /**< Dependency graph. irG[k] is the row of
                     *   G[k]. */
    const int *jcG; /**< Dependency graph. Index to data of first
                     *   non-zero element in row k. */
    const int *irS; /**< State-change matrix. irS[k] is the row of
                     *   S[k]. */
    const int *jcS; /**< State-change matrix. Index to data of first
                     *   non-zero element in row k. */
    const int *prS; /**< State-change matrix. Value of item (i, j)
                     *   in S. */
    const int *irE; /**< Select matrix for events. irE[k] is the row
                     *   of E[k]. */
    const int *jcE; /**< Select matrix for events. Index to data of
                     *   first non-zero element in row k. */

    /*** Callbacks ***/
    TRFun *tr_fun;  /**< Vector of function pointers to
                       n *   transition rate functions */
    PTSFun pts_fun; /**< Callback after each time step */

    /*** Keep track of time ***/
    double tt;           /**< The global time. */
    double next_day;     /**< The global time of next day. */
    const double *tspan; /**< Output times. tspan[0] is the start time
                          *   and tspan[length(tspan)-1] is the stop
                          *   time.*/
    int tlen;            /**< Number of sampling points in time. */
    int U_it;            /**< Index to next time in tspan */
    int V_it;            /**< Index to next time in tspan */

    /*** Data vectors ***/
    int *u;           /**< Vector with the number of individuals in
                       *   each compartment in each node in the
                       *   thread. */
    int *U;           /**< If the solution is written to a dense
                       *   matrix the compartment output is a matrix U
                       *   ((Nn * Nc) X length(tspan)). U(:,j)
                       *   contains the state of the system at
                       *   tspan(j). */
    const int *irU;   /**< If the solution is written to a sparse
                       *   matrix, irU[k] is the row of U[k]. */
    const int *jcU;   /**< If the solution is written to a sparse
                       *   matrix, index to data of first non-zero
                       *   element in row k. */
    double    *prU;   /**< If the solution is written to a sparse
                       *   matrix, value of item (i, j) in U. */
    double *v;        /**< Vector with the continuous state in each
                       *   node in the thread. */
    double *v_new;    /**< Vector with the continuous state in each
                       *   node in the thread after the post time step
                       *   function. */
    double *V;        /**< If the solution is written to a dense
                       *   matrix the continuous output is a matrix V
                       *   ((Nn * Nd) X length(tspan)). V(:,j)
                       *   contains the state of the system at
                       *   tspan(j). */
    const int *irV;   /**< If the solution is written to a sparse
                       *   matrix, irV[k] is the row of V[k]. */
    const int *jcV;   /**< If the solution is written to a sparse
                       *   matrix, index to data of first non-zero
                       *   element in row k. */
    double    *prV;   /**< If the solution is written to a sparse
                       *   matrix, value of item (i, j) in V. */
    const double *ldata; /**< Matrix (Nld X Nn). ldata(:,j) gives a
                          *   local data vector for node #j. */
    const double *gdata; /**< The global data vector. */
    const int *N;     /**< Shift matrix for internal and external
                       *   transfer events. */
    int *update_node; /**< Vector of length Nn used to indicate nodes
                       *   for update. */

    double *sum_t_rate; /**< Vector of length Nn with the sum of
                         *   propensities in every node. */
    double *t_rate;     /**< Transition rate matrix (Nt X Nn) with all
                         *   propensities for state transitions. */
    double *t_time;     /**< Time for next event (transition) in each
                         *   node. */
    int errcode;        /**< The error state of the thread. 0 if
                         *   ok. */
    gsl_rng *rng;       /**< The random number generator. */

    /*** Scheduled events ***/
    SimInf_scheduled_events *E1; /**< E1 events to process. */
    int E1_index;         /**< Index to the next E1 event to
                           *   process. */
    SimInf_scheduled_events *E2; /**< E2 events to process. */
    int E2_index;         /**< Index to the next E2 event to
                           *   process. */

    /*** Vectors for sampling individuals ***/
    int *individuals;     /**< Vector to store the result of the
                           *   sampling during scheduled events
                           *   processing. */
    int *u_tmp;           /**< Temporary vector with the compartment
                           *   state in a node when sampling
                           *   individuals for scheduled events. */

    /* (START) AEM specific variables */
    gsl_rng **rng_vec;       /**< The random number generator. */

    int *reactHeap;      /**< The reaction heap (AEM solver). */
    int reactHeapSize;

    int *reactNode;

    double *reactTimes, *reactInf; /**< holds the times (AEM) */
    /* (END) AEM specific variables */

} SimInf_thread_args;

/* Shared variables */
int n_thread_aem = 0;
int *uu_aem = NULL;
double *vv_1_aem = NULL;
double *vv_2_aem = NULL;
int *update_node_aem = NULL;
SimInf_thread_args *sim_args_aem = NULL;

/**
 * Allocate memory for scheduled events
 *
 * @param e SimInf_scheduled_events structure for events.
 * @param n Number of events.
 * @return 0 on success else SIMINF_ERR_ALLOC_MEMORY_BUFFER
 */
static int SimInf_allocate_events(SimInf_scheduled_events *e, int n)
{
    if (e && n > 0) {
        e->len = n;
        e->event = malloc(n * sizeof(int));
        if (!e->event)
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        e->time = malloc(n * sizeof(int));
        if (!e->time)
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        e->node = malloc(n * sizeof(int));
        if (!e->node)
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        e->dest = malloc(n * sizeof(int));
        if (!e->dest)
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        e->n = malloc(n * sizeof(int));
        if (!e->n)
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        e->proportion = malloc(n * sizeof(double));
        if (!e->proportion)
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        e->select = malloc(n * sizeof(int));
        if (!e->select)
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        e->shift = malloc(n * sizeof(int));
        if (!e->shift)
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
    }

    return 0;
}

/**
 * Free allocated memory to scheduled events
 *
 * @param e The SimInf_scheduled_events events to free.
 */
static void SimInf_free_events(SimInf_scheduled_events *e)
{
    if (e) {
        if (e->event)
            free(e->event);
        e->event = NULL;
        if (e->time)
            free(e->time);
        e->time = NULL;
        if (e->node)
            free(e->node);
        e->node = NULL;
        if (e->dest)
            free(e->dest);
        e->dest = NULL;
        if (e->n)
            free(e->n);
        e->n = NULL;
        if (e->proportion)
            free(e->proportion);
        e->proportion = NULL;
        if (e->select)
            free(e->select);
        e->select = NULL;
        if (e->shift)
            free(e->shift);
        e->shift = NULL;
        free(e);
    }
}

/**
 * Free allocated memory to siminf thread arguments
 */
static void SimInf_free_args(SimInf_thread_args *sa)
{
    if (sa) {
        if (sa->rng)
            gsl_rng_free(sa->rng);
        sa->rng = NULL;
        if (sa->t_rate)
            free(sa->t_rate);
        sa->t_rate = NULL;
        if (sa->sum_t_rate)
            free(sa->sum_t_rate);
        sa->sum_t_rate = NULL;
        if (sa->t_time)
            free(sa->t_time);
        sa->t_time = NULL;
        if (sa->individuals)
            free(sa->individuals);
        sa->individuals = NULL;
        if (sa->u_tmp)
            free(sa->u_tmp);
        sa->u_tmp = NULL;
        if (sa->E1)
            SimInf_free_events(sa->E1);
        sa->E1 = NULL;
        if (sa->E2)
            SimInf_free_events(sa->E2);
        sa->E2 = NULL;
        /* AEM variables */
	if(sa->rng_vec){
            for(int i = 0; i < (sa->Nn)*(sa->Nt); i++)
                gsl_rng_free(sa->rng_vec[i]);
        }
        sa->rng_vec = NULL;
	if(sa->reactHeap)
            free(sa->reactHeap);
        sa->reactHeap = NULL;
	if(sa->reactInf)
            free(sa->reactInf);
        sa->reactInf = NULL;
	if(sa->reactNode)
            free(sa->reactNode);
        sa->reactNode = NULL;
	if(sa->reactTimes)
            free(sa->reactTimes);
        sa->reactTimes = NULL;
    }
}



/**
 * Split scheduled events to E1 and E2 events by number of threads
 * used during simulation
 *
 * Thread id 0 is the main thread. All E2 events are assigned to
 * thread id 0.
 *
 * All E1 events for a node are assigned to the same thread.
 *
 * @param len Number of scheduled events.
 * @param event The type of event i.
 * @param time The time of event i.
 * @param node The source node index (one based) of event i.
 * @param dest The dest node index (one-based) of event i.
 * @param n The number of individuals in event i. n[i] >= 0.
 * @param proportion If n[i] equals zero, then the number of
 *        individuals to sample is calculated by summing the number of
 *        individuals in the states determined by select[i] and
 *        multiplying with the proportion. 0 <= p[i] <= 1.
 * @param select Column j (one-based) in the event matrix that
 *        determines the states to sample from.
 * @param shift Column j (one-based) in the shift matrix S that
 *        determines the shift of the internal and external
 *        transfer event.
 * @param Nn Total number of nodes.
 * @return 0 if Ok, else error code.
 */
static int SimInf_split_events(
    int len, const int *event, const int *time, const int *node,
    const int *dest, const int *n, const double *proportion,
    const int *select, const int *shift, int Nn)
{
    int i;
    int errcode = 0;
    int chunk_size = Nn / n_thread_aem;
    int *E1_i = NULL;
    int E2_i = 0;

    /* Split events to each thread */
    E1_i = calloc(n_thread_aem, sizeof(int));
    if (!E1_i) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    for (i = 0; i < len; i++) {
        int k;

        switch (event[i]) {
        case EXIT_EVENT:
        case ENTER_EVENT:
        case INTERNAL_TRANSFER_EVENT:
            k = (node[i] - 1) / chunk_size;
            if (k >= n_thread_aem)
                k = n_thread_aem - 1;
            E1_i[k]++;
            break;
        case EXTERNAL_TRANSFER_EVENT:
            E2_i++;
            break;
        default:
            errcode = SIMINF_UNDEFINED_EVENT;
            goto cleanup;
        }
    }

    /* Allocate memory for E1 and E2 events. */
    for (i = 0; i < n_thread_aem; i++) {
        errcode = SimInf_allocate_events(sim_args_aem[i].E1, E1_i[i]);
        if (errcode)
            goto cleanup;
        E1_i[i] = 0;

        if (i == 0) {
            errcode = SimInf_allocate_events(sim_args_aem[0].E2, E2_i);
            if (errcode)
                goto cleanup;
            E2_i = 0;
        }
    }

    for (i = 0; i < len; i++) {
        int j, k;
        SimInf_scheduled_events *e;

        switch (event[i]) {
        case EXIT_EVENT:
        case ENTER_EVENT:
        case INTERNAL_TRANSFER_EVENT:
            k = (node[i] - 1) / chunk_size;
            if (k >= n_thread_aem)
                k = n_thread_aem - 1;
            j = E1_i[k]++;
            e = sim_args_aem[k].E1;
            break;
        case EXTERNAL_TRANSFER_EVENT:
            j = E2_i++;
            e = sim_args_aem[0].E2;
            break;
        default:
            errcode = SIMINF_UNDEFINED_EVENT;
            goto cleanup;
        }

        e->event[j]      = event[i];
        e->time[j]       = time[i];
        e->node[j]       = node[i] - 1;
        e->dest[j]       = dest[i] - 1;
        e->n[j]          = n[i];
        e->proportion[j] = proportion[i];
        e->select[j]     = select[i] - 1;
        e->shift[j]      = shift[i] - 1;
    }

cleanup:
    if (E1_i)
        free(E1_i);

    return errcode;
}

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
static int SimInf_solver()
{
    int k;

#pragma omp parallel
    {
        int i;

#pragma omp for
        for (i = 0; i < n_thread_aem; i++) {
            int node;
            SimInf_thread_args sa = *&sim_args_aem[i];

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

	    *&sim_args_aem[i] = sa;
        }
    }

    /* Check for error during initialization. */
    for (k = 0; k < n_thread_aem; k++)
        if (sim_args_aem[k].errcode)
            return sim_args_aem[k].errcode;

    /* Main loop. */
    for (;;) {
#pragma omp parallel
        {
            int i;

#pragma omp for
            for (i = 0; i < n_thread_aem; i++) {
                int node;
                SimInf_thread_args sa = *&sim_args_aem[i];
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
                            uu_aem[e1.node[j] * sa.Nc + sa.irE[sa.jcE[s]]] += e1.n[j];
                            if (uu_aem[e1.node[j] * sa.Nc + sa.irE[sa.jcE[s]]] < 0)
                                sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                        }
                    } else {
                        sa.errcode = SimInf_sample_select(
                            sa.irE, sa.jcE, sa.Nc, uu_aem, e1.node[j],
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
                                uu_aem[kk] -= sa.individuals[jj];
                                if (uu_aem[kk] < 0) {
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
                                uu_aem[kk + ll] += sa.individuals[jj];
                                if (uu_aem[kk + ll] < 0) {
                                    sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                                    break;
                                }

                                /* Remove individuals from previous
                                 * compartments in node */
                                uu_aem[kk] -= sa.individuals[jj];
                                if (uu_aem[kk] < 0) {
                                    sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                                    break;
                                }
                            }
                        }
                    }

                    /* Indicate node for update */
                    update_node_aem[e1.node[j]] = 1;
                    sa.E1_index++;
                }

                *&sim_args_aem[i] = sa;
	    }


#pragma omp barrier

#pragma omp master
            {
                SimInf_thread_args sa = *&sim_args_aem[0];
                SimInf_scheduled_events e2 = *sa.E2;

                /* (3) Incorporate all scheduled E2 events */
                while (sa.E2_index < e2.len &&
                       sa.tt >= e2.time[sa.E2_index] &&
                       !sa.errcode)
                {
                    sa.errcode = SimInf_sample_select(
                        sa.irE, sa.jcE, sa.Nc, uu_aem, e2.node[sa.E2_index],
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
                        uu_aem[k1 + ll] += sa.individuals[jj];
                        if (uu_aem[k1 + ll] < 0) {
                            sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                            break;
                        }

                        /* Remove individuals from node */
                        uu_aem[k2] -= sa.individuals[jj];
                        if (uu_aem[k2] < 0) {
                            sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                            break;
                        }
                    }

                    /* Indicate node and dest for update */
                    update_node_aem[e2.node[sa.E2_index]] = 1;
                    update_node_aem[e2.dest[sa.E2_index]] = 1;
                    sa.E2_index++;
                }

                *&sim_args_aem[0] = sa;
            }

#pragma omp barrier

#pragma omp for
            for (i = 0; i < n_thread_aem; i++) {
                int node;
                SimInf_thread_args sa = *&sim_args_aem[i];

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

                *&sim_args_aem[i] = sa;
            }
        }

        /* 6b) Handle the case where the solution is stored in a sparse
         * matrix */
        while (!sim_args_aem[0].U && sim_args_aem[0].U_it < sim_args_aem[0].tlen &&
               sim_args_aem[0].tt > sim_args_aem[0].tspan[sim_args_aem[0].U_it]) {
            int j;

            /* Copy compartment state to U_sparse */
            for (j = sim_args_aem[0].jcU[sim_args_aem[0].U_it];
                 j < sim_args_aem[0].jcU[sim_args_aem[0].U_it + 1]; j++)
                sim_args_aem[0].prU[j] = sim_args_aem[0].u[sim_args_aem[0].irU[j]];
            sim_args_aem[0].U_it++;
        }

        while (!sim_args_aem[0].V && sim_args_aem[0].V_it < sim_args_aem[0].tlen &&
               sim_args_aem[0].tt > sim_args_aem[0].tspan[sim_args_aem[0].V_it]) {
            int j;

            /* Copy continuous state to V_sparse */
            for (j = sim_args_aem[0].jcV[sim_args_aem[0].V_it];
                 j < sim_args_aem[0].jcV[sim_args_aem[0].V_it + 1]; j++)
                sim_args_aem[0].prV[j] = sim_args_aem[0].v_new[sim_args_aem[0].irV[j]];
            sim_args_aem[0].V_it++;
        }

        /* Swap the pointers to the continuous state variable so that
         * 'v' equals 'v_new'. Moreover, check for error. */
        for (k = 0; k < n_thread_aem; k++) {
            double *v_tmp = sim_args_aem[k].v;
            sim_args_aem[k].v = sim_args_aem[k].v_new;
            sim_args_aem[k].v_new = v_tmp;
            if (sim_args_aem[k].errcode)
                return sim_args_aem[k].errcode;
        }

        /* If the simulation has reached the final time, exit. */
        if (sim_args_aem[0].U_it >= sim_args_aem[0].tlen)
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

    n_thread_aem = args->Nthread;

    /* Set compartment state to the initial state. */
    uu_aem = malloc(args->Nn * args->Nc * sizeof(int));
    if (!uu_aem) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    memcpy(uu_aem, args->u0, args->Nn * args->Nc * sizeof(int));

    /* Copy u0 to either U[, 1] or U_sparse[, 1] */
    if (args->U) {
        memcpy(args->U, args->u0, args->Nn * args->Nc * sizeof(int));
    } else {
        for (i = args->jcU[0]; i < args->jcU[1]; i++)
            args->prU[i] = args->u0[args->irU[i]];
    }

    /* Set continuous state to the initial state in each node. */
    vv_1_aem = malloc(args->Nn * args->Nd * sizeof(double));
    vv_2_aem = malloc(args->Nn * args->Nd * sizeof(double));
    if (!vv_1_aem || !vv_2_aem) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    memcpy(vv_1_aem, args->v0, args->Nn * args->Nd * sizeof(double));

    /* Copy v0 to either V[, 1] or V_sparse[, 1] */
    if (args->V) {
        memcpy(args->V, args->v0, args->Nn * args->Nd * sizeof(double));
    } else {
        for (i = args->jcV[0]; i < args->jcV[1]; i++)
            args->prV[i] = args->v0[args->irV[i]];
    }

    /* Setup vector to keep track of nodes that must be updated due to
     * scheduled events */
    update_node_aem = calloc(args->Nn, sizeof(int));
    if (!update_node_aem) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!rng) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    gsl_rng_set(rng, args->seed);

    sim_args_aem = calloc(n_thread_aem, sizeof(SimInf_thread_args));
    if (!sim_args_aem) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }


    for (i = 0; i < n_thread_aem; i++) {
        /* Random number generator */
        sim_args_aem[i].rng = gsl_rng_alloc(gsl_rng_mt19937);
        if (!sim_args_aem[i].rng) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }
        gsl_rng_set(sim_args_aem[i].rng, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));

        /* Constants */
        sim_args_aem[i].Ntot = args->Nn;
        sim_args_aem[i].Ni = i * (args->Nn / n_thread_aem);
        sim_args_aem[i].Nn = args->Nn / n_thread_aem;
        if (i == (n_thread_aem - 1))
            sim_args_aem[i].Nn += (args->Nn % n_thread_aem);
        sim_args_aem[i].Nt = args->Nt;
        sim_args_aem[i].Nc = args->Nc;
        sim_args_aem[i].Nd = args->Nd;
        sim_args_aem[i].Nld = args->Nld;

        /* Sparse matrices */
        sim_args_aem[i].irG = args->irG;
        sim_args_aem[i].jcG = args->jcG;
        sim_args_aem[i].irS = args->irS;
        sim_args_aem[i].jcS = args->jcS;
        sim_args_aem[i].prS = args->prS;
        sim_args_aem[i].irE = args->irE;
        sim_args_aem[i].jcE = args->jcE;

        /* Callbacks */
        sim_args_aem[i].tr_fun = args->tr_fun;
        sim_args_aem[i].pts_fun = args->pts_fun;

        /* Keep track of time */
        sim_args_aem[i].tt = args->tspan[0];
        sim_args_aem[i].next_day = floor(sim_args_aem[i].tt) + 1.0;
        sim_args_aem[i].tspan = args->tspan;
        sim_args_aem[i].tlen = args->tlen;
        sim_args_aem[i].U_it = 1;
        sim_args_aem[i].V_it = 1;

        /* start AEM specific
         */
        /* Binary heap storing all reaction events */
        /* we have one for each node. Heap is thus only the size of the # transitions */
        sim_args_aem[i].reactHeapSize = sim_args_aem[i].Nt;

        sim_args_aem[i].reactNode = (int *)malloc(sim_args_aem[i].Nn*sim_args_aem[i].Nt*sizeof(int));
        sim_args_aem[i].reactHeap = (int *)malloc(sim_args_aem[i].Nn*sim_args_aem[i].Nt*sizeof(int));
        sim_args_aem[i].reactTimes = (double *)malloc(sim_args_aem[i].Nn*sim_args_aem[i].Nt*sizeof(double));
        /* Initialize with zero elements */
        sim_args_aem[i].reactInf = (double *)calloc(sim_args_aem[i].Nn*sim_args_aem[i].Nt,sizeof(double));

        /* random generator for sample select with 1 per transition in each node */
        sim_args_aem[i].rng_vec = (gsl_rng **)malloc(sim_args_aem[i].Nn*sim_args_aem[i].Nt*sizeof(gsl_rng*));

        for(int node = 0; node < sim_args_aem[i].Nn; node++){
            /* sim_args_aem[i].rng_vec[node] = (gsl_rng **)malloc(sim_args_aem[i].Nt*sizeof(gsl_rng*)); */
            for(int trans = 0; trans < sim_args_aem[i].Nt; trans++){
                /* Random number generator */
                sim_args_aem[i].rng_vec[sim_args_aem[i].Nt*node + trans] = gsl_rng_alloc(gsl_rng_mt19937);
                if (!sim_args_aem[i].rng_vec[sim_args_aem[i].Nt*node + trans]) {
                    errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
                    goto cleanup;
                }
                gsl_rng_set(sim_args_aem[i].rng_vec[sim_args_aem[i].Nt*node + trans],
                            gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
            }
        }
        /* end AEM specific
         */

        /* Data vectors */
        sim_args_aem[i].N = args->N;
        if (args->U) {
            sim_args_aem[i].U = args->U;
        } else if (i == 0) {
            sim_args_aem[i].irU = args->irU;
            sim_args_aem[i].jcU = args->jcU;
            sim_args_aem[i].prU = args->prU;
        }
        sim_args_aem[i].u = &uu_aem[sim_args_aem[i].Ni * args->Nc];
        if (args->V) {
            sim_args_aem[i].V = args->V;
        } else if (i == 0) {
            sim_args_aem[i].irV = args->irV;
            sim_args_aem[i].jcV = args->jcV;
            sim_args_aem[i].prV = args->prV;
        }
        sim_args_aem[i].v = &vv_1_aem[sim_args_aem[i].Ni * args->Nd];
        sim_args_aem[i].v_new = &vv_2_aem[sim_args_aem[i].Ni * args->Nd];
        sim_args_aem[i].ldata = &(args->ldata[sim_args_aem[i].Ni * sim_args_aem[i].Nld]);
        sim_args_aem[i].gdata = args->gdata;
        sim_args_aem[i].update_node = &update_node_aem[sim_args_aem[i].Ni];

        /* Scheduled events */
        sim_args_aem[i].E1 = calloc(1, sizeof(SimInf_scheduled_events));
        if (!sim_args_aem[i].E1) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        if (i == 0) {
            sim_args_aem[i].E2 = calloc(1, sizeof(SimInf_scheduled_events));
            if (!sim_args_aem[i].E2) {
                errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
                goto cleanup;
            }
        }

        sim_args_aem[i].individuals = calloc(args->Nc, sizeof(int));
        if (!sim_args_aem[i].individuals) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        sim_args_aem[i].u_tmp = calloc(args->Nc, sizeof(int));
        if (!sim_args_aem[i].u_tmp) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

	/* Create transition rate matrix (Nt X Nn) and total rate
         * vector. In t_rate we store all propensities for state
         * transitions, and in sum_t_rate the sum of propensities
         * in every node. */
        sim_args_aem[i].t_rate = malloc(args->Nt * sim_args_aem[i].Nn * sizeof(double));
        if (!sim_args_aem[i].t_rate) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }
        sim_args_aem[i].sum_t_rate = malloc(sim_args_aem[i].Nn * sizeof(double));
        if (!sim_args_aem[i].sum_t_rate) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }
        sim_args_aem[i].t_time = malloc(sim_args_aem[i].Nn * sizeof(double));
        if (!sim_args_aem[i].t_time) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }
    }

    /* Split scheduled events into E1 and E2 events. */
    errcode = SimInf_split_events(
        args->len, args->event, args->time, args->node, args->dest, args->n,
        args->proportion, args->select, args->shift, args->Nn);
    if (errcode)
        goto cleanup;

    errcode = SimInf_solver();

cleanup:
    if (uu_aem) {
        free(uu_aem);
        uu_aem = NULL;
    }

    if (vv_1_aem) {
        free(vv_1_aem);
        vv_1_aem = NULL;
    }

    if (vv_2_aem) {
        free(vv_2_aem);
        vv_2_aem = NULL;
    }

    if (update_node_aem) {
        free(update_node_aem);
        update_node_aem = NULL;
    }

    if (rng)
        gsl_rng_free(rng);

    if (sim_args_aem) {
        for (i = 0; i < n_thread_aem; i++)
            SimInf_free_args(&sim_args_aem[i]);
        free(sim_args_aem);
        sim_args_aem = NULL;
    }

    return errcode;
}
