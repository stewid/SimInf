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

#include "SimInf.h"
#include "SimInf_solver.h"

/**
 * Sample individuals from a node
 *
 * Individuals are sampled from the states determined by select.
 *
 * @param irE Select matrix for events. irE[k] is the row of E[k].
 * @param jcE Select matrix for events. Index to data of first
 *        non-zero element in row k.
 * @param Nc Number of compartments in each node.
 * @param u The state vector with number of individuals in each
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
 * @param u_tmp Help vector for sampling individuals.
 * @param rng Random number generator.
 * @return 0 if Ok, else error code.
 */
static int SimInf_sample_select(
    const int *irE, const int *jcE, int Nc, const int *u,
    int node, int select, int n, double proportion,
    int *individuals, int *u_tmp, gsl_rng *rng)
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
    if (Nstates <= 0 ||     /* No states to sample from, we shouldn't be here. */
        n > Nindividuals || /* Can not sample this number of individuals       */
        n < 0)              /* Can not sample negative number of individuals.  */
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
        /* Randomly sample n individuals from Nindividudals in
         * the Nstates */
        memcpy(u_tmp, &u[node * Nc], Nc * sizeof(int));
        while (n > 0) {
            double cum, rand = gsl_rng_uniform_pos(rng) * Nindividuals;

            /* Determine from which compartment the individual was
             * sampled from */
            for (i = jcE[select], cum = u_tmp[irE[i]];
                 i < jcE[select + 1] && rand > cum;
                 i++, cum += u_tmp[irE[i]]);

            /* Update sampled individual */
            u_tmp[irE[i]]--;
            individuals[irE[i]]++;

            Nindividuals--;
            n--;
        }
    }

    return 0;
}

/**
 * Allocate memory for scheduled events
 *
 * @param e scheduled_events structure for events.
 * @param n Number of events.
 * @return 0 on success else SIMINF_ERR_ALLOC_MEMORY_BUFFER
 */
static int SimInf_allocate_events(SimInf_scheduled_events_data *e, int n)
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
 * @param Nthread Number of threads to use during simulation.
 * @return 0 if Ok, else error code.
 */
static int SimInf_split_events(
    SimInf_scheduled_events *out,
    int len, const int *event, const int *time, const int *node,
    const int *dest, const int *n, const double *proportion,
    const int *select, const int *shift, int Nn, int Nthread)
{
    int i;
    int error = 0;
    int chunk_size = Nn / Nthread;
    int *E1_i = NULL;
    int E2_i = 0;

    /* Split events to each thread */
    E1_i = calloc(Nthread, sizeof(int));
    if (!E1_i) {
        error = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
        goto cleanup;
    }

    for (i = 0; i < len; i++) {
        int k;

        switch (event[i]) {
        case EXIT_EVENT:
        case ENTER_EVENT:
        case INTERNAL_TRANSFER_EVENT:
            k = (node[i] - 1) / chunk_size;
            if (k >= Nthread)
                k = Nthread - 1;
            E1_i[k]++;
            break;
        case EXTERNAL_TRANSFER_EVENT:
            E2_i++;
            break;
        default:
            error = SIMINF_UNDEFINED_EVENT;
            goto cleanup;
        }
    }

    /* Allocate memory for E1 and E2 events. */
    for (i = 0; i < Nthread; i++) {
        error = SimInf_allocate_events(out[i].E1, E1_i[i]);
        if (error)
            goto cleanup;
        E1_i[i] = 0;

        if (i == 0) {
            error = SimInf_allocate_events(out[0].E2, E2_i);
            if (error)
                goto cleanup;
            E2_i = 0;
        }
    }

    for (i = 0; i < len; i++) {
        int j, k;
        SimInf_scheduled_events_data *e;

        switch (event[i]) {
        case EXIT_EVENT:
        case ENTER_EVENT:
        case INTERNAL_TRANSFER_EVENT:
            k = (node[i] - 1) / chunk_size;
            if (k >= Nthread)
                k = Nthread - 1;
            j = E1_i[k]++;
            e = out[k].E1;
            break;
        case EXTERNAL_TRANSFER_EVENT:
            j = E2_i++;
            e = out[0].E2;
            break;
        default:
            error = SIMINF_UNDEFINED_EVENT;
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

    return error;
}

/**
 * Create and initialize data to process scheduled events. The
 * generated data structure must be freed by the user.
 *
 * @param out the resulting data structure.
 * @param args structure with data for the solver.
 * @param rng random number generator
 * @return 0 or an error code
 */
int SimInf_scheduled_events_create(
    SimInf_scheduled_events **out, SimInf_solver_args *args, gsl_rng *rng)
{
    int error = SIMINF_ERR_ALLOC_MEMORY_BUFFER, i;
    SimInf_scheduled_events *events = NULL;

    events = calloc(args->Nthread, sizeof(SimInf_scheduled_events));
    if (!events)
        goto on_error;

    for (i = 0; i < args->Nthread; i++) {
        /* Matrices to process events */
        events[i].irE = args->irE;
        events[i].jcE = args->jcE;
        events[i].N = args->N;

        /* Scheduled events */
        events[i].E1 = calloc(1, sizeof(SimInf_scheduled_events_data));
        if (!events[i].E1)
            goto on_error;

        if (i == 0) {
            events[i].E2 = calloc(1, sizeof(SimInf_scheduled_events_data));
            if (!events[i].E2)
                goto on_error;
        }

        events[i].individuals = calloc(args->Nc, sizeof(int));
        if (!events[i].individuals)
            goto on_error;

        events[i].u_tmp = calloc(args->Nc, sizeof(int));
        if (!events[i].u_tmp)
            goto on_error;

        /* Random number generator */
        events[i].rng = gsl_rng_alloc(gsl_rng_mt19937);
        if (!events[i].rng)
            goto on_error;
        gsl_rng_set(events[i].rng, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
    }

    /* Split scheduled events into E1 and E2 events. */
    error = SimInf_split_events(
        events, args->len, args->event, args->time, args->node,
        args->dest, args->n, args->proportion, args->select,
        args->shift, args->Nn, args->Nthread);
    if (error)
        goto on_error;

    *out = events;
    return 0;

on_error:
    SimInf_scheduled_events_free(events, args->Nthread);
    return error;
}

/**
 * Free allocated memory for scheduled events data
 *
 * @param e SimInf_scheduled_events_data to free
 */
static void SimInf_scheduled_events_data_free(
    SimInf_scheduled_events_data *events)
{
    if (events) {
        if (events->event)
            free(events->event);
        events->event = NULL;
        if (events->time)
            free(events->time);
        events->time = NULL;
        if (events->node)
            free(events->node);
        events->node = NULL;
        if (events->dest)
            free(events->dest);
        events->dest = NULL;
        if (events->n)
            free(events->n);
        events->n = NULL;
        if (events->proportion)
            free(events->proportion);
        events->proportion = NULL;
        if (events->select)
            free(events->select);
        events->select = NULL;
        if (events->shift)
            free(events->shift);
        events->shift = NULL;
        free(events);
    }
}

/**
 * Free allocated memory to process events
 *
 * @param events SimInf_scheduled_events to free
 * @param Nthread number of threads that was used during simulation.
 */
void SimInf_scheduled_events_free(
    SimInf_scheduled_events *events, int Nthread)
{
    if (events) {
        int i;

        for (i = 0; i < Nthread; i++) {
            SimInf_scheduled_events *e = &events[i];

            if (e) {
                if (e->E1)
                    SimInf_scheduled_events_data_free(e->E1);
                e->E1 = NULL;
                if (e->E2)
                    SimInf_scheduled_events_data_free(e->E2);
                e->E2 = NULL;
                if (e->individuals)
                    free(e->individuals);
                e->individuals = NULL;
                if (e->u_tmp)
                    free(e->u_tmp);
                e->u_tmp = NULL;
                if (e->rng)
                    gsl_rng_free(e->rng);
                e->rng = NULL;
            }
        }

        free(events);
    }
}

void SimInf_process_E1_events(
    SimInf_compartment_model *model,
    SimInf_scheduled_events *events)
{
    SimInf_compartment_model m = *&model[0];
    SimInf_scheduled_events e = *&events[0];
    SimInf_scheduled_events_data e1 = *e.E1;

    while (e.E1_index < e1.len &&
           m.tt >= e1.time[e.E1_index] &&
           !m.error)
    {
        const int i = e.E1_index;
        const int s = e1.select[i];
        const int node = e1.node[i] - m.Ni;

        if (e1.node[i] < 0 || e1.node[i] >= m.Ntot) {
            m.error = SIMINF_ERR_NODE_OUT_OF_BOUNDS;
            break;
        }

        if (e1.event[i] == ENTER_EVENT) {
            /* All individuals enter first non-zero compartment,
             * i.e. a non-zero entry in element in the select
             * column. */
            if (e.jcE[s] < e.jcE[s + 1]) {
                m.u[node * m.Nc + e.irE[e.jcE[s]]] += e1.n[i];
                if (m.u[node * m.Nc + e.irE[e.jcE[s]]] < 0)
                    m.error = SIMINF_ERR_NEGATIVE_STATE;
            }
        } else {
            m.error = SimInf_sample_select(
                e.irE, e.jcE, m.Nc, m.u, node,
                e1.select[i], e1.n[i], e1.proportion[i],
                e.individuals, e.u_tmp, e.rng);

            if (m.error)
                break;

            if (e1.event[i] == EXIT_EVENT) {
                int ii;

                for (ii = e.jcE[s]; ii < e.jcE[s + 1]; ii++) {
                    const int jj = e.irE[ii];
                    const int kk = node * m.Nc + jj;

                    /* Remove individuals from node */
                    m.u[kk] -= e.individuals[jj];
                    if (m.u[kk] < 0) {
                        m.error = SIMINF_ERR_NEGATIVE_STATE;
                        break;
                    }
                }
            } else { /* INTERNAL_TRANSFER_EVENT */
                int ii;

                for (ii = e.jcE[s]; ii < e.jcE[s + 1]; ii++) {
                    const int jj = e.irE[ii];
                    const int kk = node * m.Nc + jj;
                    const int ll = e.N[e1.shift[i] * m.Nc + jj];

                    /* Add individuals to new compartments in node */
                    m.u[kk + ll] += e.individuals[jj];
                    if (m.u[kk + ll] < 0) {
                        m.error = SIMINF_ERR_NEGATIVE_STATE;
                        break;
                    }

                    /* Remove individuals from previous compartments
                     * in node */
                    m.u[kk] -= e.individuals[jj];
                    if (m.u[kk] < 0) {
                        m.error = SIMINF_ERR_NEGATIVE_STATE;
                        break;
                    }
                }
            }
        }

        /* Indicate node for update */
        m.update_node[node] = 1;
        e.E1_index++;
    }

    *&events[0] = e;
    *&model[0] = m;
}

void SimInf_process_E2_events(
    SimInf_compartment_model *model,
    SimInf_scheduled_events *events)
{
    SimInf_compartment_model m = *&model[0];
    SimInf_scheduled_events e = *&events[0];
    SimInf_scheduled_events_data e2 = *e.E2;

    /* Incorporate all scheduled E2 events */
    while (e.E2_index < e2.len &&
           m.tt >= e2.time[e.E2_index] &&
           !m.error)
    {
        int i;
        const int dest = e2.dest[e.E2_index];
        const int node = e2.node[e.E2_index];

        if (dest < 0 || dest >= m.Ntot) {
            m.error = SIMINF_ERR_DEST_OUT_OF_BOUNDS;
            break;
        }

        if (node < 0 || node >= m.Ntot) {
            m.error = SIMINF_ERR_NODE_OUT_OF_BOUNDS;
            break;
        }

        m.error = SimInf_sample_select(
            e.irE, e.jcE, m.Nc, m.u, node,
            e2.select[e.E2_index], e2.n[e.E2_index],
            e2.proportion[e.E2_index], e.individuals,
            e.u_tmp, e.rng);

        if (m.error)
            break;

        for (i = e.jcE[e2.select[e.E2_index]];
             i < e.jcE[e2.select[e.E2_index] + 1];
             i++)
        {
            const int jj = e.irE[i];
            const int k1 = dest * m.Nc + jj;
            const int k2 = node * m.Nc + jj;
            const int ll = e2.shift[e.E2_index] < 0 ? 0 :
                e.N[e2.shift[e.E2_index] * m.Nc + jj];

            /* Add individuals to dest */
            m.u[k1 + ll] += e.individuals[jj];
            if (m.u[k1 + ll] < 0) {
                m.error = SIMINF_ERR_NEGATIVE_STATE;
                break;
            }

            /* Remove individuals from node */
            m.u[k2] -= e.individuals[jj];
            if (m.u[k2] < 0) {
                m.error = SIMINF_ERR_NEGATIVE_STATE;
                break;
            }
        }

        /* Indicate node and dest for update */
        m.update_node[node] = 1;
        m.update_node[dest] = 1;
        e.E2_index++;
    }

    *&events[0] = e;
    *&model[0] = m;
}

/**
 * Handle the case where the solution is stored in a sparse matrix
 *
 * Store solution if tt has passed the next time in tspan. Report
 * solution up to, but not including tt.
 *
 * @param SimInf_compartment_model *model data to store.
 */
void SimInf_store_solution_sparse(SimInf_compartment_model *model)
{
    while (!model[0].U && model[0].U_it < model[0].tlen &&
           model[0].tt > model[0].tspan[model[0].U_it]) {
        int j;

        /* Copy compartment state to U_sparse */
        for (j = model[0].jcU[model[0].U_it];
             j < model[0].jcU[model[0].U_it + 1]; j++)
            model[0].prU[j] = model[0].u[model[0].irU[j]];
        model[0].U_it++;
    }

    while (!model[0].V && model[0].V_it < model[0].tlen &&
           model[0].tt > model[0].tspan[model[0].V_it]) {
        int j;

        /* Copy continuous state to V_sparse */
        for (j = model[0].jcV[model[0].V_it];
             j < model[0].jcV[model[0].V_it + 1]; j++)
            model[0].prV[j] = model[0].v_new[model[0].irV[j]];
        model[0].V_it++;
    }
}

/**
 * Free allocated memory for an epidemiological compartment
 * model.
 *
 * @param model the data structure to free.
 * @param Nthread number of threads that was used during simulation.
 */
void SimInf_compartment_model_free(SimInf_compartment_model *model, int Nthread)
{
    if (model) {
        int i;

        for (i = 0; i < Nthread; i++) {
            SimInf_compartment_model *m = &model[i];

            if (m) {
                free(m->t_rate);
                m->t_rate = NULL;
                free(m->sum_t_rate);
                m->sum_t_rate = NULL;
                free(m->t_time);
                m->t_time = NULL;
            }
        }

        free(model[0].u);
        model[0].u = NULL;
        free(model[0].v);
        model[0].v = NULL;
        free(model[0].v_new);
        model[0].v_new = NULL;
        free(model[0].update_node);
        model[0].update_node = NULL;
        free(model);
    }
}

/**
 * Create and initialize data for an epidemiological compartment
 * model. The generated model must be freed by the user.
 *
 * @param out the resulting data structure.
 * @param args structure with data for the solver.
 * @return 0 or SIMINF_ERR_ALLOC_MEMORY_BUFFER
 */
int SimInf_compartment_model_create(
    SimInf_compartment_model **out, SimInf_solver_args *args)
{
    int i;
    SimInf_compartment_model *model = NULL;

    /* Allocate memory for the compartment model. */
    model = calloc(args->Nthread, sizeof(SimInf_compartment_model));
    if (!model)
        goto on_error;

    /* Allocate memory to keep track of the continuous state in each
     * node. */
    model[0].v = malloc(args->Nn * args->Nd * sizeof(double));
    if (!model[0].v)
        goto on_error;
    model[0].v_new = malloc(args->Nn * args->Nd * sizeof(double));
    if (!model[0].v_new)
        goto on_error;

    /* Set continuous state to the initial state in each node. */
    memcpy(model[0].v, args->v0, args->Nn * args->Nd * sizeof(double));

    /* Setup vector to keep track of nodes that must be updated due to
     * scheduled events */
    model[0].update_node = calloc(args->Nn, sizeof(int));
    if (!model[0].update_node)
        goto on_error;

    /* Allocate memory for compartment state and set compartment state
     * to the initial state. */
    model[0].u = malloc(args->Nn * args->Nc * sizeof(int));
    if (!model[0].u)
        goto on_error;
    memcpy(model[0].u, args->u0, args->Nn * args->Nc * sizeof(int));

    for (i = 0; i < args->Nthread; i++) {
        /* Constants */
        model[i].Ntot = args->Nn;
        model[i].Ni = i * (args->Nn / args->Nthread);
        model[i].Nn = args->Nn / args->Nthread;
        if (i == (args->Nthread - 1))
            model[i].Nn += (args->Nn % args->Nthread);
        model[i].Nt = args->Nt;
        model[i].Nc = args->Nc;
        model[i].Nd = args->Nd;
        model[i].Nld = args->Nld;

        /* Sparse matrices */
        model[i].irG = args->irG;
        model[i].jcG = args->jcG;
        model[i].irS = args->irS;
        model[i].jcS = args->jcS;
        model[i].prS = args->prS;

        /* Callbacks */
        model[i].tr_fun = args->tr_fun;
        model[i].pts_fun = args->pts_fun;

        /* Keep track of time */
        model[i].tt = args->tspan[0];
        model[i].next_unit_of_time = floor(model[i].tt) + 1.0;
        model[i].tspan = args->tspan;
        model[i].tlen = args->tlen;
        model[i].U_it = 1;
        model[i].V_it = 1;

        /* Data vectors */
        if (args->U) {
            model[i].U = args->U;
        } else if (i == 0) {
            model[i].irU = args->irU;
            model[i].jcU = args->jcU;
            model[i].prU = args->prU;
        }

        if (args->V) {
            model[i].V = args->V;
        } else if (i == 0) {
            model[i].irV = args->irV;
            model[i].jcV = args->jcV;
            model[i].prV = args->prV;
        }

        if (i > 0) {
            model[i].u = &model[0].u[model[i].Ni * args->Nc];
            model[i].v = &model[0].v[model[i].Ni * args->Nd];
            model[i].v_new = &model[0].v_new[model[i].Ni * args->Nd];
            model[i].update_node = &model[0].update_node[model[i].Ni];
        }

        model[i].ldata = &(args->ldata[model[i].Ni * model[i].Nld]);
        model[i].gdata = args->gdata;

        /* Create transition rate matrix (Nt X Nn) and total rate
         * vector. In t_rate we store all propensities for state
         * transitions, and in sum_t_rate the sum of propensities in
         * every node. */
        model[i].t_rate = malloc(args->Nt * model[i].Nn * sizeof(double));
        if (!model[i].t_rate)
            goto on_error;
        model[i].sum_t_rate = malloc(model[i].Nn * sizeof(double));
        if (!model[i].sum_t_rate)
            goto on_error;
        model[i].t_time = malloc(model[i].Nn * sizeof(double));
        if (!model[i].t_time)
            goto on_error;
    }

    *out = model;

    return 0;

on_error:
    SimInf_compartment_model_free(model, args->Nthread);
    return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
}
