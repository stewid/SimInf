/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015 - 2016 Stefan Engblom
 *  Copyright (C) 2015 - 2016 Stefan Widgren
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

#include "siminf.h"

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
 * Structure that represents scheduled events.
 */
typedef struct scheduled_events
{
    int len;            /**< Number of events. */
    int *event;         /**< The type of event i. */
    int *time;          /**< The time of event i. */
    int *node;          /**< The source node of event i. */
    int *dest;          /**< The dest node of event i. */
    int *n;             /**< The number of individuals in the scheduled
                         *   event. n[i] >= 0. */
    double *proportion; /**< If n[i] equals zero, then the number of
                         *   individuals to sample is calculated by
                         *   summing the number of individuals in the
                         *   states determined by select[i] and
                         *   multiplying with the proportion.
                         *   0 <= p[i] <= 1. */
    int *select;        /**< Column j in the event matrix that
                         *   determines the states to sample from. */
    int *shift;         /**< Column j in the shift matrix that
                         *   determines the shift of the internal
                         *   and external transfer event. */
} scheduled_events;

/* Maximum number of individuals to sample from */
#define MAX_INDIVIDUALS 10000

/**
 * Structure to hold thread specific data/arguments for simulation.
 */
typedef struct siminf_thread_args
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
                     *   transition rate functions */
    PTSFun pts_fun; /**< Callback after each time step */

    /*** Keep track of time ***/
    double tt;           /**< The global time. */
    double next_day;     /**< The global time of next day. */
    const double *tspan; /**< Output times. tspan[0] is the start time
                          *   and tspan[length(tspan)-1] is the stop
                          *   time.*/
    int tlen;            /**< Number of sampling points in time. */
    int it;              /**< Index to next time in tspan */

    /*** Data vectors ***/
    int *U;           /**< The compartment output is a matrix U ((Nn *
                       *   Nc) X length(tspan)). U(:,j) contains the
                       *   state of the system at tspan(j). */
    int *u;           /**< Vector with compartment state in each node
                       *   in thread. */
    double *V;        /**< The continuous output is a matrix V ((Nn *
                       *   Nd) X length(tspan)). V(:,j) contains the
                       *   state of the system at tspan(j). */
    double *v;        /**< Vector with continuous state in each node
                       *   in thread. */
    double *v_new;    /**< Vector with continuous state in each node
                       *   in thread after post time step function. */
    const double *ldata; /**< Matrix (Nld X Nn). ldata(:,j) gives a
                          *   local data vector for node #j. */
    const double *gdata; /**< The global data vector. */
    const int *sd;    /**< Each node can be assigned to a
                       *   sub-domain. */
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
    scheduled_events *E1; /**< E1 events to process. */
    int E1_index;         /**< Index to the next E1 event to
                           *   process. */
    scheduled_events *E2; /**< E2 events to process. */
    int E2_index;         /**< Index to the next E2 event to
                           *   process. */

    /*** Vectors for sampling individuals ***/
    int *individuals;               /**< Vector to store the result of
                                     *   the sampling during scheduled
                                     *   events processing. */
    int kind[MAX_INDIVIDUALS];      /**< Help vector for sampling
                                     *   individuals. */
    int kind_dest[MAX_INDIVIDUALS]; /**< Help vector for sampling
                                     *   individuals. */
} siminf_thread_args;

/* Shared variables */
int n_thread = 0;
int *uu = NULL;
double *vv_1 = NULL;
double *vv_2 = NULL;
int *update_node = NULL;
siminf_thread_args *sim_args = NULL;

/**
 * Allocate memory for scheduled events
 *
 * @param e scheduled_events structure for events.
 * @param n Number of events.
 * @return 0 on success else SIMINF_ERR_ALLOC_MEMORY_BUFFER
 */
static int siminf_allocate_events(scheduled_events *e, int n)
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
 * @param e The scheduled_events events to free.
 */
static void siminf_free_events(scheduled_events *e)
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
static void siminf_free_args(siminf_thread_args *sa)
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
        if (sa->E1)
            siminf_free_events(sa->E1);
        sa->E1 = NULL;
        if (sa->E2)
            siminf_free_events(sa->E2);
        sa->E2 = NULL;
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
static int siminf_split_events(
    int len, const int *event, const int *time, const int *node,
    const int *dest, const int *n, const double *proportion,
    const int *select, const int *shift, int Nn)
{
    int i;
    int errcode = 0;
    int chunk_size = Nn / n_thread;
    int *E1_i = NULL;
    int E2_i = 0;

    /* Split events to each thread */
    E1_i = calloc(n_thread, sizeof(int));
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
            if (k >= n_thread)
                k = n_thread - 1;
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
    for (i = 0; i < n_thread; i++) {
        errcode = siminf_allocate_events(sim_args[i].E1, E1_i[i]);
        if (errcode)
            goto cleanup;
        E1_i[i] = 0;

        if (i == 0) {
            errcode = siminf_allocate_events(sim_args[0].E2, E2_i);
            if (errcode)
                goto cleanup;
            E2_i = 0;
        }
    }

    for (i = 0; i < len; i++) {
        int j, k;
        scheduled_events *e;

        switch (event[i]) {
        case EXIT_EVENT:
        case ENTER_EVENT:
        case INTERNAL_TRANSFER_EVENT:
            k = (node[i] - 1) / chunk_size;
            if (k >= n_thread)
                k = n_thread - 1;
            j = E1_i[k]++;
            e = sim_args[k].E1;
            break;
        case EXTERNAL_TRANSFER_EVENT:
            j = E2_i++;
            e = sim_args[0].E2;
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
 * Sample individuals from a node
 *
 * Individuals are sampled from the states determined by select.
 *
 * @param irE Select matrix for events. irE[k] is the row of E[k].
 * @param jcE Select matrix for events. Index to data of first
 *        non-zero element in row k.
 * @param Nc Number of compartments in each node.
 * @param state The state vector with number of individuals in each
 *        compartment at each node. The current state in each node is
 *        offset by node * Nc.
 * @param node The node to sample.
 * @param select Column j in the Select matrix that determines the
 *        states to sample from.
 * @param n The number of individuals to sample. n >= 0.
 * @param proportion If n equals zero, then the number of individuals
 *        to sample is calculated by summing the number of individuals
 *        in the states determined by select and multiplying with the
 *        proportion. 0 <= proportion <= 1.
 * @param individuals The result of the sampling is stored in the
 *        individuals vector.
 * @param kind Help vector for sampling individuals.
 * @param kind_des Help vector for sampling individuals.
 * @param rng Random number generator.
 * @return 0 if Ok, else error code.
 */
static int sample_select(
    const int *irE, const int *jcE, int Nc, const int *u,
    int node, int select, int n, double proportion,
    int *individuals, int *kind, int *kind_dest, gsl_rng *rng)
{
    int i, Nstates, Nindividuals = 0, Nkinds = 0;

    /* Clear vector with number of sampled individuals */
    memset(individuals, 0, Nc * sizeof(int));

    /* 1) Count number of states with individuals */
    /* 2) Count total number of individuals       */
    for (i = jcE[select]; i < jcE[select + 1]; i++) {
        int nk = u[node * Nc + irE[i]];
        if (nk > 0)
            Nkinds++;
        Nindividuals += nk;
    }

    /* Number of states */
    Nstates = jcE[select + 1] - jcE[select];

    /* If n == 0, use the proportion of Nindividuals, else use n as */
    /* the number of individuals to sample                          */
    if (n == 0)
        n = round(proportion * Nindividuals);

    /* Error checking. */
    if (Nstates <= 0        /* No states to sample from, we shouldn't be here. */
        || n > Nindividuals /* Can not sample this number of individuals       */
        || n < 0)           /* Can not sample negative number of individuals.  */
        return SIMINF_ERR_SAMPLE_SELECT;

    /* Handle cases that require no random sampling */
    if (n == 0) {
        /* We are done */
        return 0;
    } else if (Nindividuals == n) {
        /* Include all individuals */
        for (i = jcE[select]; i < jcE[select + 1]; i++)
            individuals[irE[i]] = u[node * Nc + irE[i]];
        return 0;
    } else if (Nstates == 1) {
        /* Only individuals from one state to select from. */
        individuals[irE[jcE[select]]] = n;
        return 0;
    } else if (Nkinds == 1) {
        /* All individuals to choose from in one state */
        for (i = jcE[select]; i < jcE[select + 1]; i++) {
            if (u[node * Nc + irE[i]] > 0) {
                individuals[irE[i]] = n;
                break;
            }
        }
        return 0;
    }

    /* Handle cases that require random sampling */
    if (Nstates == 2) {
        /* Sample from the hypergeometric distribution */
        i = jcE[select];
        individuals[irE[i]] = gsl_ran_hypergeometric(
            rng,
            u[node * Nc + irE[i]],
            u[node * Nc + irE[i+1]],
            n);
        individuals[irE[i+1]] = n - individuals[irE[i]];
    } else {
        /* Randomly choose n individuals from a vector of
         * Nindividudals in Nstates */
        int j;

        /* Intialize and populate kind vector */
        if (Nindividuals > MAX_INDIVIDUALS)
            return 1;
        for (i = jcE[select], j = 0; i < jcE[select + 1]; i++) {
            int k, nk, l;

            k  = irE[i];           /* The kind  */
            nk = u[node * Nc + k]; /* N of kind */

            /* Set kind 'k' for 'nk' individuals */
            for (l = 0; l < nk; l++)
                kind[j++] = k;
        }

        /* Randomly choose n individuals from kind vector */
        gsl_ran_choose(rng, kind_dest, n, kind, Nindividuals, sizeof(int));

        /* Count kind of the choosen individuals */
        for (i = 0; i < n; i++)
            individuals[kind_dest[i]]++;
    }

    return 0;
}

/**
 * Siminf solver
 *
 * @return 0 if Ok, else error code.
 */
static int siminf_solver()
{
    #pragma omp parallel
    {
        int i;

        #pragma omp for
        for (i = 0; i < n_thread; i++) {
            int node;
            siminf_thread_args sa = *&sim_args[i];

            /* Initialize transition rate and time to event. Calculate
             * the transition rate for every transition and every
             * node. Store the sum of the transition rates in each
             * node in sum_t_rate. Calculate time to next event
             * (transition) in each node. */
            for (node = 0; node < sa.Nn; node++) {
                int j;

                sa.sum_t_rate[node] = 0.0;
                for (j = 0; j < sa.Nt; j++) {
                    sa.t_rate[node * sa.Nt + j] =
                        (*sa.tr_fun[j])(&sa.u[node * sa.Nc],
                                        &sa.v[node * sa.Nd],
                                        &sa.ldata[node * sa.Nld],
                                        sa.gdata,
                                        sa.tt,
                                        sa.sd[node]);

                    sa.sum_t_rate[node] += sa.t_rate[node * sa.Nt + j];
                }

                sa.t_time[node] =
                    -log(1.0 - gsl_rng_uniform(sa.rng)) /
                    sa.sum_t_rate[node] + sa.tt;
            }

            *&sim_args[i] = sa;
        }
    }

    /* Main loop. */
    for (;;) {
        int k;

        #pragma omp parallel
        {
            int i;

            #pragma omp for
            for (i = 0; i < n_thread; i++) {
                int node;
                siminf_thread_args sa = *&sim_args[i];
                scheduled_events e1 = *sa.E1;

                /* (1) Handle internal epidemiological model,
                 * continuous-time Markov chain. */
                for (node = 0; node < sa.Nn && !sa.errcode; node++) {
                    while (sa.t_time[node] < sa.next_day) {
                        double cum, rand, tot_rate, delta = 0.0;
                        int j, tr;

                        /* a) Determine the transition that did occur
                         * (directSSA). */
                        rand = gsl_rng_uniform_pos(sa.rng) * sa.sum_t_rate[node];
                        for (tr = 0, cum = sa.t_rate[node * sa.Nt];
                             tr < sa.Nt && rand > cum;
                             tr++, cum += sa.t_rate[node * sa.Nt + tr]);

                        /* b) Update the state of the node */
                        for (j = sa.jcS[tr]; j < sa.jcS[tr + 1]; j++) {
                            sa.u[node * sa.Nc + sa.irS[j]] += sa.prS[j];
                            if (sa.u[node * sa.Nc + sa.irS[j]] < 0)
                                sa.errcode = SIMINF_ERR_NEGATIVE_STATE;
                        }

                        /* c) Recalculate sum_t_rate[node] using
                         * dependency graph. */
                        for (j = sa.jcG[tr]; j < sa.jcG[tr + 1]; j++) {
                            const double old = sa.t_rate[node * sa.Nt + sa.irG[j]];
                            delta += (sa.t_rate[node * sa.Nt + sa.irG[j]] =
                                      (*sa.tr_fun[sa.irG[j]])(
                                          &sa.u[node * sa.Nc],
                                          &sa.v[node * sa.Nd],
                                          &sa.ldata[node * sa.Nld],
                                          sa.gdata,
                                          sa.t_time[node],
                                          sa.sd[node])) - old;
                        }
                        sa.sum_t_rate[node] += delta;

                        /* d) Compute time to new event for this node. */
                        tot_rate = sa.sum_t_rate[node];
                        if (tot_rate > 0.0) {
                            sa.t_time[node] =
                                -log(1.0 - gsl_rng_uniform(sa.rng)) /
                                tot_rate + sa.t_time[node];
                        } else {
                            sa.t_time[node] = INFINITY;
                        }
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
                        sa.errcode = sample_select(
                            sa.irE, sa.jcE, sa.Nc, uu, e1.node[j],
                            e1.select[j], e1.n[j], e1.proportion[j],
                            sa.individuals, sa.kind, sa.kind_dest, sa.rng);

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
                siminf_thread_args sa = *&sim_args[0];
                scheduled_events e2 = *sa.E2;

                /* (3) Incorporate all scheduled E2 events */
                while (sa.E2_index < e2.len &&
                       sa.tt >= e2.time[sa.E2_index] &&
                       !sa.errcode)
                {
                    sa.errcode = sample_select(
                        sa.irE, sa.jcE, sa.Nc, uu, e2.node[sa.E2_index],
                        e2.select[sa.E2_index], e2.n[sa.E2_index],
                        e2.proportion[sa.E2_index], sa.individuals,
                        sa.kind, sa.kind_dest, sa.rng);

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
            for (i = 0; i < n_thread; i++) {
                int node;
                siminf_thread_args sa = *&sim_args[i];

                /* (4) Incorporate model specific actions after each
                 * timestep e.g. update the infectious pressure
                 * variable. Moreover, update transition rates in
                 * nodes that are indicated for update */
                for (node = 0; node < sa.Nn; node++) {
                    const int rc = sa.pts_fun(
                        &sa.v_new[node * sa.Nd], &sa.u[node * sa.Nc],
                        &sa.v[node * sa.Nd], &sa.ldata[node * sa.Nld],
                        sa.gdata, sa.Ni + node, sa.tt, sa.sd[node]);

                    if (rc < 0) {
                        sa.errcode = rc;
                        break;
                    } else if (rc > 0 || sa.update_node[node]) {
                        /* Update transition rates */
                        int j = 0;
                        double delta = 0.0, old_t_rate = sa.sum_t_rate[node];

                        for (; j < sa.Nt; j++) {
                            const double old = sa.t_rate[node * sa.Nt + j];
                            delta += (sa.t_rate[node * sa.Nt + j] =
                                      (*sa.tr_fun[j])(
                                          &sa.u[node * sa.Nc],
                                          &sa.v_new[node * sa.Nd],
                                          &sa.ldata[node * sa.Nld],
                                          sa.gdata, sa.tt, sa.sd[node])) - old;
                        }
                        sa.sum_t_rate[node] += delta;

                        if (sa.sum_t_rate[node] > 0.0) {
                            if (old_t_rate > 0.0 && !isinf(old_t_rate)) {
                                sa.t_time[node] =
                                    old_t_rate / sa.sum_t_rate[node]
                                    * (sa.t_time[node] - sa.tt) + sa.tt;
                            } else {
                                sa.t_time[node] =
                                    -log(1.0 - gsl_rng_uniform(sa.rng)) /
                                    sa.sum_t_rate[node] + sa.tt;
                            }
                        } else {
                            sa.t_time[node] = INFINITY;
                        }

                        sa.update_node[node] = 0;
                    }
                }

                /* (5) The global time now equals next_day. */
                sa.tt = sa.next_day;
                sa.next_day += 1.0;

                /* (6) Store solution if tt has passed the next time
                 * in tspan. Report solution up to, but not including
                 * tt. */
                if (sa.tt > sa.tspan[sa.it]) {
                    while (sa.it < sa.tlen && sa.tt > sa.tspan[sa.it]) {
                        /* Copy compartment state to U */
                        memcpy(&sa.U[sa.Nc * (sa.Ntot * sa.it + sa.Ni)],
                               sa.u, sa.Nn * sa.Nc * sizeof(int));

                        /* Copy continuous state to V */
                        memcpy(&sa.V[sa.Nd * ((sa.Ntot * sa.it++) + sa.Ni)],
                               sa.v_new, sa.Nn * sa.Nd * sizeof(double));
                    }
                }

                *&sim_args[i] = sa;
            }
        }

        /* Swap the pointers to the continuous state variable so that
         * 'v' equals 'v_new'. Moreover, check for error. */
        for (k = 0; k < n_thread; k++) {
            double *v_tmp = sim_args[k].v;
            sim_args[k].v = sim_args[k].v_new;
            sim_args[k].v_new = v_tmp;
            if (sim_args[k].errcode)
                return sim_args[k].errcode;
        }

        /* If the simulation has reached the final time, exit. */
        if (sim_args[0].it >= sim_args[0].tlen)
            break;
    }

    return 0;
}

/**
 * Initialize and run siminf solver
 *
 * G is a sparse matrix dependency graph (Nt X Nt) in compressed
 * column format (CCS). A non-zeros entry in element i of column j
 * indicates that transition rate i needs to be recalculated if the
 * transition j occurs.
 *
 * S is the state-changing sparse matrix (Nc X Nt) in compressed
 * column format (CCS). Each column corresponds to a transition, and
 * execution of transition j amounts to adding the j'th column to the
 * state vector.
 *
 * @param u0 Initial state vector u0. Integer (Nc X Nn). Gives the
 *        initial number of individuals in each compartment in every
 *        node.
 * @param v0 Initial continuous state vector v0. Double (Nd X Nn).
 *        Gives the initial value of the continuous state variables
 *        in every node.
 * @param irG Dependency graph. irG[k] is the row of G[k].
 * @param jcG Dependency graph. Index to data of first non-zero
 *        element in row k.
 * @param irS State-change matrix. irS[k] is the row of S[k].
 * @param jcS State-change matrix. Index to data of first non-zero
 *        element in row k.
 * @param prS State-change matrix. Value of item (i, j) in S.
 * @param tspan Double vector. Output times. tspan[0] is the start
 *        time and tspan[length(tspan)-1] is the stop time.
 * @param tlen Number of sampling points in time.
 * @param U The output is a matrix U ((Nn * Nc) X length(tspan)).
 *        U(:,j) contains the state of the system at tspan(j).
 * @param V The continuous state output is a matrix V
 *        ((Nn * Nd) X length(tspan)).
 *        V(:,j) contains the continuous state of the system at tspan(j).
 * @param ldata Double matrix (Nld X Nn). Generalized data matrix,
 *        data(:,j) gives a local data vector for node #j.
 * @param gdata The global data vector.
 * @param sd Integer vector of length Nn. Each node can be assigned to
 *        a sub-domain.
 * @param Nn Number of nodes.
 * @param Nc Number of compartments in each node.
 * @param Nt Total number of different transitions.
 * @param Nd Number of continuous state variables.
 * @param Nld Length of the local data vector 'ldata' for each
          node. The 'ldata' vector is sent to propensities and the
          post time step function.
 * @param irE Select matrix for events. irE[k] is the row of E[k].
 * @param jcE Select matrix for events. Index to data of first
 *        non-zero element in row k.
 * @param N Shift matrix for internal transfer events.
 * @param len Number of events.
 * @param event The type of event i.
 * @param time The time of event i.
 * @param node The source node of event i.
 * @param dest The dest node of event i.
 * @param n The number of individuals in the scheduled event. n[i] >= 0.
 * @param proportion If n[i] equals zero, then the number of
 *        individuals to sample is calculated by summing the number of
 *        individuals in the states determined by select[i] and
 *        multiplying with the proportion. 0 <= p[i] <= 1.
 * @param select Column j in the event matrix E that determines the
 *        states to sample from.
 * @param shift Column j in the shift matrix S that determines the
 *        shift of the internal and external transfer event.
 * @param Nthread Number of threads to use during simulation.
 * @param seed Random number seed.
 * @param tr_fun Vector of function pointers to transition rate functions.
 * @param pts_fun Function pointer to callback after each time step
 *        e.g. to update the infectious pressure.
 * @return 0 if Ok, else error code.
 */
int siminf_run_solver(
    const int *u0, const double *v0, const int *irG, const int *jcG,
    const int *irS, const int *jcS, const int *prS, const double *tspan,
    int tlen, int *U, double *V, const double *ldata, const double *gdata,
    const int *sd, int Nn, int Nc, int Nt, int Nd, int Nld, const int *irE,
    const int *jcE, const int *N, int len, const int *event,
    const int *time, const int *node, const int *dest, const int *n,
    const double *proportion, const int *select, const int *shift,
    int Nthread, unsigned long int seed, TRFun *tr_fun, PTSFun pts_fun)
{
    int i, errcode;
    gsl_rng *rng = NULL;

#ifdef _OPENMP
    if (Nthread < 1)
        Nthread = omp_get_num_procs();
#else
    Nthread = 1;
#endif
    if (Nn < Nthread)
        n_thread = Nn;
    else
        n_thread = Nthread;
#ifdef _OPENMP
    omp_set_num_threads(n_thread);
#endif

    /* Set compartment state to the initial state. */
    uu = malloc(Nn * Nc * sizeof(int));
    if (!uu) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    memcpy(uu, u0, Nn * Nc * sizeof(int));

    /* Set continuous state to the initial state in each node. */
    vv_1 = malloc(Nn * Nd * sizeof(double));
    vv_2 = malloc(Nn * Nd * sizeof(double));
    if (!vv_1 || !vv_2) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    memcpy(vv_1, v0, Nn * Nd * sizeof(double));

    /* Setup vector to keep track of nodes that must be updated due to
     * scheduled events */
    update_node = calloc(Nn, sizeof(int));
    if (!update_node) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!rng) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }
    gsl_rng_set(rng, seed);

    sim_args = calloc(n_thread, sizeof(siminf_thread_args));
    if (!sim_args) {
        errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    for (i = 0; i < n_thread; i++) {
        /* Random number generator */
        sim_args[i].rng = gsl_rng_alloc(gsl_rng_mt19937);
        if (!sim_args[i].rng) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }
        gsl_rng_set(sim_args[i].rng, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));

        /* Constants */
        sim_args[i].Ntot = Nn;
        sim_args[i].Ni = i * (Nn / n_thread);
        sim_args[i].Nn = Nn / n_thread;
        if (i == (n_thread - 1))
            sim_args[i].Nn += (Nn % n_thread);
        sim_args[i].Nt = Nt;
        sim_args[i].Nc = Nc;
        sim_args[i].Nd = Nd;
        sim_args[i].Nld = Nld;

        /* Sparse matrices */
        sim_args[i].irG = irG;
        sim_args[i].jcG = jcG;
        sim_args[i].irS = irS;
        sim_args[i].jcS = jcS;
        sim_args[i].prS = prS;
        sim_args[i].irE = irE;
        sim_args[i].jcE = jcE;

        /* Callbacks */
        sim_args[i].tr_fun = tr_fun;
        sim_args[i].pts_fun = pts_fun;

        /* Keep track of time */
        sim_args[i].tt = tspan[0];
        sim_args[i].next_day = floor(sim_args[i].tt) + 1.0;
        sim_args[i].tspan = tspan;
        sim_args[i].tlen = tlen;

        /* Data vectors */
        sim_args[i].N = N;
        sim_args[i].U = U;
        sim_args[i].u = &uu[sim_args[i].Ni * Nc];
        sim_args[i].V = V;
        sim_args[i].v = &vv_1[sim_args[i].Ni * Nd];
        sim_args[i].v_new = &vv_2[sim_args[i].Ni * Nd];
        sim_args[i].ldata = &ldata[sim_args[i].Ni * sim_args[i].Nld];
        sim_args[i].gdata = gdata;
        sim_args[i].sd = &sd[sim_args[i].Ni];
        sim_args[i].update_node = &update_node[sim_args[i].Ni];

        /* Scheduled events */
        sim_args[i].E1 = calloc(1, sizeof(scheduled_events));
        if (!sim_args[i].E1) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        if (i == 0) {
            sim_args[i].E2 = calloc(1, sizeof(scheduled_events));
            if (!sim_args[i].E2) {
                errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
                goto cleanup;
            }
        }

        sim_args[i].individuals = calloc(Nc, sizeof(int));
        if (!sim_args[i].individuals) {
            errcode = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
            goto cleanup;
        }

        /* Create transition rate matrix (Nt X Nn) and total rate
         * vector. In t_rate we store all propensities for state
         * transitions, and in sum_t_rate the sum of propensities
         * in every node. */
        sim_args[i].t_rate = malloc(Nt * sim_args[i].Nn * sizeof(double));
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
    errcode = siminf_split_events(
        len, event, time, node, dest, n, proportion, select, shift, Nn);
    if (errcode)
        goto cleanup;

    errcode = siminf_solver();

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
        for (i = 0; i < n_thread; i++)
            siminf_free_args(&sim_args[i]);
        free(sim_args);
        sim_args = NULL;
    }

    return errcode;
}
