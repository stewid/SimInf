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

#ifdef _OPENMP
#include <omp.h>
#endif

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
 * @param prE Select matrix for events. Value of item E[i, j].
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
 *        in the states determined by select and sampling the from a
 *        binomial distribution using the proportion and the number
 *        of individuals in the compartments. 0 <= proportion <= 1.
 * @param individuals The result of the sampling is stored in the
 *        individuals vector.
 * @param rng Random number generator.
 * @return 0 if Ok, else error code.
 */
static int
SimInf_sample_select(
    const int *irE,
    const int *jcE,
    const double *prE,
    int Nc,
    const int *u,
    int node,
    int select,
    int n,
    double proportion,
    int *individuals,
    gsl_rng *rng)
{
    int i, Nstates, Nindividuals = 0, Nkinds = 0;

    /* Clear vector with number of sampled individuals */
    memset(individuals, 0, Nc * sizeof(int));

    /* 1) Count number of states with individuals */
    /* 2) Count total number of individuals       */
    for (i = jcE[select]; i < jcE[select + 1]; i++) {
        const int nk = u[node * Nc + irE[i]];
        if (nk > 0)
            Nkinds++;
        Nindividuals += nk;
    }

    /* Number of states */
    Nstates = jcE[select + 1] - jcE[select];

    /* If n == 0, use the proportion of Nindividuals, else use n as */
    /* the number of individuals to sample.                         */
    if (n == 0) {
        if (proportion < 0 || proportion > 1)
            return SIMINF_ERR_INVALID_PROPORTION;
        n = (int)gsl_ran_binomial(rng, proportion, Nindividuals);
    }

    /* Error checking. */
    if (Nstates <= 0 ||     /* No states to sample from, we shouldn't be here. */
        n > Nindividuals || /* Cannot sample this number of individuals.       */
        n < 0)              /* Cannot sample negative number of individuals.   */
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

    /* Determine if all weights are identical in the column in E and
     * the individuals can be sampled from a hypergeometric
     * distribution or if they need to be sampled from a biased
     * urn. */
    for (i = jcE[select] + 1; i < jcE[select + 1]; i++) {
        if (prE[i] != prE[i - 1])
            goto sample_biased_urn;
    }

    /* All weights are equal. Sample from the hypergeometric
     * distribution. For a multivariate hypergeometric distribution,
     * use the algortihm described by James E. Gentle (2003, page 206)
     * in 'Random Number Generation and Monte Carlo Methods'.*/
    for (i = jcE[select]; i < jcE[select + 1] - 1; i++) {
        if (n == 0)
            break;

        individuals[irE[i]] = (int)gsl_ran_hypergeometric(
            rng, u[node * Nc + irE[i]],
            Nindividuals - u[node * Nc + irE[i]], n);

        Nindividuals -= u[node * Nc + irE[i]];
        n -= individuals[irE[i]];
    }

    individuals[irE[i]] = n;

    return 0;

sample_biased_urn:
    /* Non-equal weights. Perform sampling from Wallenius' noncentral
     * hypergeometric distribution by simulating an urn experiment
     * with bias and without replacement. The probability of taking an
     * individual from a compartment at a particular draw is equal to
     * this compartment’s fraction of the total weight of all
     * individuals that lie in the urn at this moment. The
     * implementation below is based on 'Simulating the urn
     * Experiment' in Fog (2008) 'Sampling Methods for Wallenius' and
     * Fisher's Noncentral Hypergeometric
     * Distributions'. Communications In statictics, Simulation and
     * Computation, 2008, vol. 37, no. 2, pp. 241-257. */

    /* Repeat the sampling until all n individuals have been taken. */
    while (n > 0) {
        double rand, cum = 0;

        /* Determine the total weight. */
        for (i = jcE[select]; i < jcE[select + 1]; i++)
            cum += prE[i] * (u[node * Nc + irE[i]] - individuals[irE[i]]);

        /* Use inversion to determine the compartment that was
         * sampled. */
        rand = gsl_rng_uniform_pos(rng) * cum;
        for (i = jcE[select], cum = prE[i] * (u[node * Nc + irE[i]] - individuals[irE[i]]);
             i < jcE[select + 1] && rand > cum;
             i++, cum += prE[i] * (u[node * Nc + irE[i]] - individuals[irE[i]]));

        /* Elaborate floating point fix: */
        if (i >= jcE[select + 1])
            i = jcE[select + 1] - 1;
        if ((prE[i] * (u[node * Nc + irE[i]] - individuals[irE[i]])) == 0.0) {
            /* Go backwards and try to find the first nonzero
             * compartment */
            for (;
                 i > jcE[select] && (prE[i] * (u[node * Nc + irE[i]] - individuals[irE[i]])) == 0.0;
                 i--);

            if ((prE[i] * (u[node * Nc + irE[i]] - individuals[irE[i]])) == 0.0) {
                /* No nonzero compartment found. */
                return SIMINF_ERR_SAMPLE_SELECT;
            }
        }

        /* Add the sampled individual. */
        individuals[irE[i]] += 1;
        n--;
    }

    return 0;
}

/**
 * Sample individuals to enter a node
 *
 * Individuals are sampled from the states determined by select.
 *
 * @param irE Select matrix for events. irE[k] is the row of E[k].
 * @param jcE Select matrix for events. Index to data of first
 *        non-zero element in row k.
 * @param prE Select matrix for events. Value of item E[i, j].
 * @param Nc Number of compartments in each node.
 * @param u The state vector with number of individuals in each
 *        compartment at each node. The current state in each node is
 *        offset by node * Nc.
 * @param node The node to sample.
 * @param select Column j in the select matrix that determines the
 *        states to enter individuals in.
 * @param n The number of individuals to enter. n >= 0.
 * @param proportion If n equals zero, then the number of individuals
 *        to enter is calculated by summing the number of individuals
 *        in the states determined by select and sampling the from a
 *        binomial distribution using the proportion and the number
 *        of individuals in the compartments. 0 <= proportion <= 1.
 * @param individuals The result of the sampling is stored in the
 *        individuals vector.
 * @param rng Random number generator.
 * @return 0 if Ok, else error code.
 */
static int
SimInf_sample_select_enter(
    const int *irE,
    const int *jcE,
    const double *prE,
    int Nc,
    const int *u,
    int node,
    int select,
    int n,
    double proportion,
    int *individuals,
    gsl_rng *rng)
{
    int i, Nstates = jcE[select + 1] - jcE[select];
    double w_cum = 0;

    /* Clear vector with number of individuals to enter */
    memset(individuals, 0, Nc * sizeof(int));

    /* If n == 0, use the proportion of individuals in the seleceted
     * compartments, else use n as the number of individuals to
     * enter. */
    if (n == 0) {
        if (proportion < 0 || proportion > 1)
            return SIMINF_ERR_INVALID_PROPORTION;

        /* Count the total number of individuals in the selected
         * compartments. */
        for (i = jcE[select]; i < jcE[select + 1]; i++)
            n += u[node * Nc + irE[i]];
        n = (int)gsl_ran_binomial(rng, proportion, n);
    }

    /* Error checking. */
    if (Nstates <= 0 || /* No compartments to enter individuals in. */
        n < 0)          /* No individuals to enter.                 */
        return SIMINF_ERR_SAMPLE_SELECT;

    /* Handle cases that require no random sampling */
    if (n == 0) {
        /* We are done */
        return 0;
    } else if (Nstates == 1) {
        /* All individuals enter one state. */
        individuals[irE[jcE[select]]] = n;
        return 0;
    }

    /* Determine the total weight. */
    for (i = jcE[select]; i < jcE[select + 1]; i++)
        w_cum += prE[i];

    /* Repeat the sampling until all n individuals have been entered. */
    while (n > 0) {
        double cum, rand = gsl_rng_uniform_pos(rng) * w_cum;

        /* Use inversion to determine the compartment that was
         * sampled. */
        for (i = jcE[select], cum = prE[i];
             i < jcE[select + 1] && rand > cum;
             i++, cum += prE[i]);

        /* Elaborate floating point fix: */
        if (i >= jcE[select + 1])
            i = jcE[select + 1] - 1;
        if (prE[i] == 0.0) {
            /* Go backwards and try to find the first nonzero
             * weight. */
            for (; i > jcE[select] && prE[i] == 0.0; i--);

            /* Check if a nonzero weight was found. */
            if (prE[i] == 0.0)
                return SIMINF_ERR_SAMPLE_SELECT;
        }

        /* Add the sampled individual. */
        individuals[irE[i]] += 1;
        n--;
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
 */
static void
SimInf_split_events(
    SimInf_scheduled_events *out,
    int len,
    const int *event,
    const int *time,
    const int *node,
    const int *dest,
    const int *n,
    const double *proportion,
    const int *select,
    const int *shift,
    int Nn,
    int Nthread)
{
    const int chunk_size = Nn / Nthread;

    for (int i = 0; i < len; i++) {
        const SimInf_scheduled_event e = {event[i], time[i], node[i] - 1,
                                          dest[i] - 1, n[i], proportion[i],
                                          select[i] - 1, shift[i] - 1};

        if (event[i] == EXTERNAL_TRANSFER_EVENT) {
            kv_push(SimInf_scheduled_event, out[0].events, e);
        } else {
            int j = (node[i] - 1) / chunk_size;
            if (j >= Nthread)
                j = Nthread - 1;
            kv_push(SimInf_scheduled_event, out[j].events, e);
        }
    }
}

/**
 * Copy scheduled events to each thread used during simulation
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
 */
static void
SimInf_copy_events(
    SimInf_scheduled_events *out,
    int len,
    const int *event,
    const int *time,
    const int *node,
    const int *dest,
    const int *n,
    const double *proportion,
    const int *select,
    const int *shift,
    int Nthread)
{
    for (int i = 0; i < len; i++) {
        const SimInf_scheduled_event e = {event[i], time[i], node[i] - 1,
                                          dest[i] - 1, n[i], proportion[i],
                                          select[i] - 1, shift[i] - 1};

        for (int j = 0; j < Nthread; j++) {
            kv_push(SimInf_scheduled_event, out[j].events, e);
        }
    }
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
attribute_hidden
int
SimInf_scheduled_events_create(
    SimInf_scheduled_events **out,
    const SimInf_solver_args *args,
    gsl_rng *rng)
{
    SimInf_scheduled_events *events = NULL;

    events = calloc(args->Nthread, sizeof(SimInf_scheduled_events));
    if (!events)
        goto on_error; /* #nocov */

    for (int i = 0; i < args->Nthread; i++) {
        /*** Constants ***/
        events[i].Nthread = args->Nthread;

        /* Matrices to process events */
        events[i].irE = args->irE;
        events[i].jcE = args->jcE;
        events[i].prE = args->prE;
        events[i].N = args->N;

        /* Scheduled events */
	kv_init(events[i].events);

        events[i].individuals = calloc(args->Nc, sizeof(int));
        if (!events[i].individuals)
            goto on_error; /* #nocov */

        /* Random number generator */
        events[i].rng = gsl_rng_alloc(gsl_rng_mt19937);
        if (!events[i].rng)
            goto on_error; /* #nocov */
        gsl_rng_set(events[i].rng, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
    }

    if (args->Nrep > 1) {
        /* Copy scheduled events into each thread. */
        SimInf_copy_events(
            events,
            args->len,
            args->event,
            args->time,
            args->node,
            args->dest,
            args->n,
            args->proportion,
            args->select,
            args->shift,
            args->Nthread);
    } else {
        /* Split scheduled events into E1 and E2 events. */
        SimInf_split_events(
            events,
            args->len,
            args->event,
            args->time,
            args->node,
            args->dest,
            args->n,
            args->proportion,
            args->select,
            args->shift,
            args->Nn,
            args->Nthread);
    }

    *out = events;
    return 0;

on_error:                                  /* #nocov */
    SimInf_scheduled_events_free(events);  /* #nocov */
    return SIMINF_ERR_ALLOC_MEMORY_BUFFER; /* #nocov */
}

/**
 * Free allocated memory to process events
 *
 * @param events SimInf_scheduled_events to free
 * @param Nthread number of threads that was used during simulation.
 */
attribute_hidden
void
SimInf_scheduled_events_free(
    SimInf_scheduled_events *events)
{
    if (events) {
        for (int i = 0; i < events->Nthread; i++) {
            SimInf_scheduled_events *e = &events[i];

            kv_destroy(e->events);
            free(e->individuals);
            e->individuals = NULL;
            gsl_rng_free(e->rng);
            e->rng = NULL;
        }

        free(events);
    }
}

/**
 * Print event information to facilitate debugging.
 *
 * @param e The SimInf_scheduled_events object to print.
 * @param irE Select matrix for events. irE[k] is the row of E[k].
 * @param jcE Select matrix for events. Index to data of first
 *        non-zero element in row k.
 * @param Nc Number of compartments in each node.
 * @param u The state vector with number of individuals in each
 *        compartment at each node. The current state in each node is
 *        offset by node * Nc.
 * @param node The node in u.
 */
static void
SimInf_print_event(
    const SimInf_scheduled_event *e,
    const int *irE,
    const int *jcE,
    const int Nc,
    const int *u,
    const int node,
    const int dest)
{
    #ifdef _OPENMP
    #  pragma omp critical
    #endif
    {
        if (irE && jcE && u) {
            int Nindividuals = 0;

            /* Count total number of individuals */
            for (int i = jcE[e->select]; i < jcE[e->select + 1]; i++)
                Nindividuals += u[node * Nc + irE[i]];

            /* Number of states */
            if ((jcE[e->select + 1] - jcE[e->select]) <= 0)
                Rprintf("No states to sample from.\n");

            if (e->n > Nindividuals)
                REprintf("Cannot sample %i for event from %i individuals.\n",
                         e->n, Nindividuals);

            if (e->n < 0)
                REprintf("Cannot sample %i individuals for event.\n", e->n);

            REprintf("\n");
        }

        if (u && (node >= 0)) {
            REprintf("Current state in node\n");
            REprintf("---------------------\n");

            REprintf("{");
            for (int i = 0; i < Nc; i++) {
                REprintf("%i", u[node * Nc + i]);
                if (i < (Nc - 1))
                    REprintf(", ");
            }
            REprintf("}\n\n");
        }

        if (u && (dest >= 0)) {
            REprintf("Current state in dest\n");
            REprintf("---------------------\n");

            REprintf("{");
            for (int i = 0; i < Nc; i++) {
                REprintf("%i", u[dest * Nc + i]);
                if (i < (Nc - 1))
                    REprintf(", ");
            }
            REprintf("}\n\n");
        }

        REprintf("Scheduled event\n");
        REprintf("---------------\n");

        switch (e->event) {
        case EXIT_EVENT:
            REprintf("event: %i (exit event)\n", e->event);
            break;
        case ENTER_EVENT:
            REprintf("event: %i (enter event)\n", e->event);
            break;
        case INTERNAL_TRANSFER_EVENT:
            REprintf("event: %i (internal transfer event)\n", e->event);
            break;
        case EXTERNAL_TRANSFER_EVENT:
            REprintf("event: %i (external transfer event)\n", e->event);
            break;
        default:
            REprintf("event: %i (undefined event)\n", e->event);
            break;
        }

        REprintf("time: %i\n", e->time);
        REprintf("node: %i\n", e->node + 1); /* One based in events data */
        REprintf("dest: %i\n", e->dest + 1); /* One based in events data */
        REprintf("n: %i\n", e->n);
        REprintf("proportion: %g\n", e->proportion);
        REprintf("select: %i\n", e->select + 1); /* One based in events data */
        REprintf("shift: %i\n\n", e->shift + 1); /* One based in events data */

        R_FlushConsole();
    }
}

/**
 * Process all scheduled E1 and E2 events where time is less or equal
 * to the global time in the simulation.
 *
 * @param model The compartment model with information for each node
 * and the global time.
 * @param events Data with events to process.
 * @param process_E2 Process only E1 events (process_E2 = 0), else
 * process both E1 and E2 events.
 * @return 0 if Ok, else error code.
 */
attribute_hidden
void
SimInf_process_events(
    SimInf_compartment_model *model,
    SimInf_scheduled_events *events,
    int process_E2)
{
    SimInf_compartment_model m = *&model[0];
    SimInf_scheduled_events e = *&events[0];

    /* Process events */
    while (e.events_index < kv_size(e.events) && !m.error) {
        const SimInf_scheduled_event ee = kv_A(e.events, e.events_index);

        if (ee.time > m.tt)
            goto done;

        if (ee.node < 0 || ee.node >= m.Ntot) {
            SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
            m.error = SIMINF_ERR_NODE_OUT_OF_BOUNDS;
            goto done;
        }

        switch (ee.event) {
        case EXIT_EVENT:
            m.error = SimInf_sample_select(
                e.irE, e.jcE, e.prE, m.Nc, m.u, ee.node - m.Ni, ee.select,
                ee.n, ee.proportion, e.individuals, e.rng);

            if (m.error) {
                SimInf_print_event(&ee, e.irE, e.jcE, m.Nc,
                                   m.u, ee.node - m.Ni, -1);
                goto done;
            }

            for (int i = e.jcE[ee.select]; i < e.jcE[ee.select + 1]; i++) {
                const int jj = e.irE[i];
                const int kn = (ee.node - m.Ni) * m.Nc + jj;

                /* Remove individuals from node */
                m.u[kn] -= e.individuals[jj];
                if (m.u[kn] < 0) {
                    SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                       m.u, ee.node - m.Ni, -1);
                    m.error = SIMINF_ERR_NEGATIVE_STATE;
                    goto done;
                }
            }
            break;

        case ENTER_EVENT:
            m.error = SimInf_sample_select_enter(
                e.irE, e.jcE, e.prE, m.Nc, m.u, ee.node - m.Ni, ee.select,
                ee.n, ee.proportion, e.individuals, e.rng);

            if (m.error) {
                SimInf_print_event(&ee, e.irE, e.jcE, m.Nc,
                                   m.u, ee.node - m.Ni, -1);
                goto done;
            }

            for (int i = e.jcE[ee.select]; i < e.jcE[ee.select + 1]; i++) {
                const int jj = e.irE[i];
                const int kn = (ee.node - m.Ni) * m.Nc + jj;

                if (ee.shift < 0) {
                    /* Add individuals to node without shifting
                     * compartments. */
                    m.u[kn] += e.individuals[jj];
                    if (m.u[kn] < 0) {
                        SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                           m.u, ee.node - m.Ni, -1);
                        m.error = SIMINF_ERR_NEGATIVE_STATE;
                        goto done;
                    }
                } else if (!e.N) {
                    /* Not possible to shift when N is not defined. */
                    SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                    m.error = SIMINF_ERR_EVENTS_N;
                    goto done;
                } else {
                    /* Process an enter event that also involves a
                     * shift between compartments. */
                    const int ll = e.N[ee.shift * m.Nc + jj];

                    /* Check that the index to the new compartment is
                     * not out of bounds. */
                    if (jj + ll < 0 || jj + ll >= m.Nc) {
                        SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                        m.error = SIMINF_ERR_SHIFT_OUT_OF_BOUNDS;
                        goto done;
                    }

                    /* Add individuals to node. */
                    m.u[kn + ll] += e.individuals[jj];
                    if (m.u[kn + ll] < 0) {
                        SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                           m.u, ee.node, -1);
                        m.error = SIMINF_ERR_NEGATIVE_STATE;
                        goto done;
                    }
                }
            }
            break;

        case INTERNAL_TRANSFER_EVENT:
            if (!e.N) {
                /* Not possible to shift when N is not defined. */
                SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                m.error = SIMINF_ERR_EVENTS_N;
                goto done;
            }

            if (ee.shift < 0) {
                /* Invalid shift parameter. */
                SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                m.error = SIMINF_ERR_EVENT_SHIFT;
                goto done;
            }

            m.error = SimInf_sample_select(
                e.irE, e.jcE, e.prE, m.Nc, m.u, ee.node - m.Ni, ee.select,
                ee.n, ee.proportion, e.individuals, e.rng);

            if (m.error) {
                SimInf_print_event(&ee, e.irE, e.jcE, m.Nc,
                                   m.u, ee.node - m.Ni, -1);
                goto done;
            }

            for (int i = e.jcE[ee.select]; i < e.jcE[ee.select + 1]; i++) {
                const int jj = e.irE[i];
                const int kn = (ee.node - m.Ni) * m.Nc + jj;
                const int ll = e.N[ee.shift * m.Nc + jj];

                /* Check that the index to the new compartment is not
                 * out of bounds. */
                if (jj + ll < 0 || jj + ll >= m.Nc) {
                    SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                    m.error = SIMINF_ERR_SHIFT_OUT_OF_BOUNDS;
                    goto done;
                }

                /* Add individuals to new compartments in node */
                m.u[kn + ll] += e.individuals[jj];
                if (m.u[kn + ll] < 0) {
                    SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                       m.u, ee.node - m.Ni, -1);
                    m.error = SIMINF_ERR_NEGATIVE_STATE;
                    goto done;
                }

                /* Remove individuals from previous compartments in
                 * node */
                m.u[kn] -= e.individuals[jj];
                if (m.u[kn] < 0) {
                    SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                       m.u, ee.node - m.Ni, -1);
                    m.error = SIMINF_ERR_NEGATIVE_STATE;
                    goto done;
                }
            }
            break;

        case EXTERNAL_TRANSFER_EVENT:
            /* Check if we are done because we only want to process E1
             * events. */
            if (!process_E2)
                goto done;

            if (ee.dest < 0 || ee.dest >= m.Ntot) {
                SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                m.error = SIMINF_ERR_DEST_OUT_OF_BOUNDS;
                goto done;
            }

            m.error = SimInf_sample_select(
                e.irE, e.jcE, e.prE, m.Nc, m.u, ee.node, ee.select, ee.n,
                ee.proportion, e.individuals, e.rng);

            if (m.error) {
                SimInf_print_event(&ee, e.irE, e.jcE, m.Nc,
                                   m.u, ee.node, ee.dest);
                goto done;
            }

            for (int i = e.jcE[ee.select]; i < e.jcE[ee.select + 1]; i++) {
                const int jj = e.irE[i];
                const int kd = ee.dest * m.Nc + jj;
                const int kn = ee.node * m.Nc + jj;

                if (ee.shift < 0) {
                    /* Add individuals to dest without shifting
                     * compartments */
                    m.u[kd] += e.individuals[jj];
                    if (m.u[kd] < 0) {
                        SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                           m.u, ee.node, ee.dest);
                        m.error = SIMINF_ERR_NEGATIVE_STATE;
                        goto done;
                    }
                } else if (!e.N) {
                    /* Not possible to shift when N is not defined. */
                    SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                    m.error = SIMINF_ERR_EVENTS_N;
                    goto done;
                } else {
                    /* Process a movement event that also involves a
                     * shift between compartments. */
                    const int ll = e.N[ee.shift * m.Nc + jj];

                    /* Check that the index to the new compartment is
                     * not out of bounds. */
                    if (jj + ll < 0 || jj + ll >= m.Nc) {
                        SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                        m.error = SIMINF_ERR_SHIFT_OUT_OF_BOUNDS;
                        goto done;
                    }

                    /* Add individuals to dest */
                    m.u[kd + ll] += e.individuals[jj];
                    if (m.u[kd + ll] < 0) {
                        SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                           m.u, ee.node, ee.dest);
                        m.error = SIMINF_ERR_NEGATIVE_STATE;
                        goto done;
                    }
                }

                /* Remove individuals from node */
                m.u[kn] -= e.individuals[jj];
                if (m.u[kn] < 0) {
                    SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                       m.u, ee.node, ee.dest);
                    m.error = SIMINF_ERR_NEGATIVE_STATE;
                    goto done;
                }
            }

            /* Indicate dest for update */
            m.update_node[ee.dest] = 1;
            break;

        default:
            SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
            m.error = SIMINF_UNDEFINED_EVENT;
            break;
        }

        /* Indicate node for update */
        m.update_node[ee.node - m.Ni] = 1;

        e.events_index++;
    }

done:
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
attribute_hidden
void
SimInf_store_solution_sparse(
    SimInf_compartment_model *model)
{
    while (!model[0].U && model[0].U_it < model[0].tlen &&
           model[0].tt > model[0].tspan[model[0].U_it]) {
        /* Copy compartment state to U_sparse */
        for (int j = model[0].jcU[model[0].U_it];
             j < model[0].jcU[model[0].U_it + 1]; j++)
            model[0].prU[j] = model[0].u[model[0].irU[j]];
        model[0].U_it++;
    }

    while (!model[0].V && model[0].V_it < model[0].tlen &&
           model[0].tt > model[0].tspan[model[0].V_it]) {
        /* Copy continuous state to V_sparse */
        for (int j = model[0].jcV[model[0].V_it];
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
attribute_hidden
void
SimInf_compartment_model_free(
    SimInf_compartment_model *model)
{
    if (model) {
        for (int i = 0; i < model->Nthread; i++) {
            SimInf_compartment_model *m = &model[i];

            free(m->t_rate);
            m->t_rate = NULL;
            free(m->sum_t_rate);
            m->sum_t_rate = NULL;
            free(m->t_time);
            m->t_time = NULL;

            if (m->Nrep > 0 || i == 0) {
                free(m->update_node);
                m->update_node = NULL;
            }
        }

        free(model[0].u);
        model[0].u = NULL;
        free(model[0].v);
        model[0].v = NULL;
        free(model[0].v_new);
        model[0].v_new = NULL;
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
attribute_hidden
int
SimInf_compartment_model_create(
    SimInf_compartment_model **out,
    SimInf_solver_args *args)
{
    SimInf_compartment_model *model = NULL;

    /* Allocate memory for the compartment model. */
    const R_xlen_t Nthread = args->Nthread;
    model = calloc(Nthread, sizeof(SimInf_compartment_model));
    if (!model)
        goto on_error; /* #nocov */

    /* Allocate memory to keep track of the continuous state in each
     * node. */
    const R_xlen_t Nrep = args->Nrep;
    const R_xlen_t Nn = args->Nn;
    const R_xlen_t Nd = args->Nd;
    model[0].v = malloc(Nrep * Nn * Nd * sizeof(double));
    if (!model[0].v)
        goto on_error; /* #nocov */
    model[0].v_new = malloc(Nrep * Nn * Nd * sizeof(double));
    if (!model[0].v_new)
        goto on_error; /* #nocov */

    /* Set continuous state to the initial state in each node. */
    memcpy(model[0].v, args->v0, Nrep * Nn * Nd * sizeof(double));
    memcpy(model[0].v_new, args->v0, Nrep * Nn * Nd * sizeof(double));

    /* Setup vector to keep track of nodes that must be updated due to
     * scheduled events */
    model[0].update_node = calloc(Nn, sizeof(int));
    if (!model[0].update_node)
        goto on_error; /* #nocov */

    /* Allocate memory for compartment state and set compartment state
     * to the initial state. */
    const R_xlen_t Nc = args->Nc;
    model[0].u = malloc(Nrep * Nn * Nc * sizeof(int));
    if (!model[0].u)
        goto on_error; /* #nocov */
    memcpy(model[0].u, args->u0, Nrep * Nn * Nc * sizeof(int));

    const R_xlen_t Nt = args->Nt;
    const R_xlen_t Nld = args->Nld;
    const R_xlen_t tlen = args->tlen;
    for (R_xlen_t i = 0; i < Nthread; i++) {
        /* Constants */
        model[i].Nthread = (int)Nthread;
        model[i].Ntot = (int)Nn;
        model[i].Nt = (int)Nt;
        model[i].Nc = (int)Nc;
        model[i].Nd = (int)Nd;
        model[i].Nld = (int)Nld;

        if (Nrep > 1) {
            /* All nodes belong to the same thread when running
             * multiple replicates of a model. */
            const R_xlen_t l = Nrep * i / Nthread;
            const R_xlen_t u = Nrep * (i + 1) / Nthread;

            model[i].Ni = 0;
            model[i].Nn = (int)Nn;
            model[i].Nrep = (int)(u - l);

            if (i > 0) {
                model[i].u = &(model[0].u[l * Nn * Nc]);
                model[i].v = &(model[0].v[l * Nn * Nd]);
                model[i].v_new = &(model[0].v_new[l * Nn * Nd]);

                /* Setup vector to keep track of nodes that must be
                 * updated due to scheduled events */
                model[i].update_node = calloc(Nn, sizeof(int));
                if (!model[i].update_node)
                    goto on_error; /* #nocov */
            }

            if (args->U)
                model[i].U = &args->U[tlen * l * Nn * Nc];
            if (args->V)
                model[i].V = &args->V[tlen * l * Nn * Nd];
        } else {
            /* The nodes are split between the threads when running
             * one replicate of a model. */
            model[i].Ni = (int)(i * (Nn / Nthread));
            model[i].Nn = (int)(Nn / Nthread);
            if (i == (Nthread - 1))
                model[i].Nn += (int)(Nn % Nthread);

            /* To ensure allocated memory in a multi-model can be
             * identified and released. */
            model[i].Nrep = 0;

            if (i > 0) {
                model[i].u = &(model[0].u[model[i].Ni * Nc]);
                model[i].v = &(model[0].v[model[i].Ni * Nd]);
                model[i].v_new = &(model[0].v_new[model[i].Ni * Nd]);
                model[i].update_node = &(model[0].update_node[model[i].Ni]);
            }

            if (args->U) {
                model[i].U = args->U;
            } else {
                model[i].irU = args->irU;
                model[i].jcU = args->jcU;
                model[i].prU = args->prU;
            }

            if (args->V) {
                model[i].V = args->V;
            } else {
                model[i].irV = args->irV;
                model[i].jcV = args->jcV;
                model[i].prV = args->prV;
            }
        }

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
        model[i].tlen = tlen;
        model[i].U_it = 0;
        model[i].V_it = 0;

        /* Data vectors */
        model[i].ldata = &(args->ldata[model[i].Ni * Nld]);
        model[i].gdata = args->gdata;

        /* Create transition rate matrix (Nt X Nn) and total rate
         * vector. In t_rate we store all propensities for state
         * transitions, and in sum_t_rate the sum of propensities in
         * every node. */
        model[i].t_rate = malloc(Nt * model[i].Nn * sizeof(double));
        if (!model[i].t_rate)
            goto on_error; /* #nocov */
        model[i].sum_t_rate = malloc(model[i].Nn * sizeof(double));
        if (!model[i].sum_t_rate)
            goto on_error; /* #nocov */
        model[i].t_time = malloc(model[i].Nn * sizeof(double));
        if (!model[i].t_time)
            goto on_error; /* #nocov */
    }

    *out = model;
    return 0;

on_error:                                  /* #nocov */
    SimInf_compartment_model_free(model);  /* #nocov */
    return SIMINF_ERR_ALLOC_MEMORY_BUFFER; /* #nocov */
}

/**
 * Print node status/information to facilitate debugging.
 *
 * @param Nc Number of compartments in each node.
 * @param u The state vector with number of individuals in each
 *        compartment in the node.
 * @param Nd Number of continuous state variables in the node.
 * @param v The continuous state vector in the node.
 * @param Nld Number of local data variables in the node.
 * @param ldata The local data vector with variables in the node.
 * @param node Zero-based index to node.
 * @param tt The current time in node.
 * @param rate The propensity. Only reported if it's infinite or less
 *        than zero.
 * @param transition Zero-based index with the state transition.
 */
attribute_hidden
void
SimInf_print_status(
    const int Nc,
    const int *u,
    const int Nd,
    const double *v,
    const int Nld,
    const double *ldata,
    const int node,
    const double tt,
    const double rate,
    const int transition)
{
    #ifdef _OPENMP
    #  pragma omp critical
    #endif
    {
        REprintf("Status:\n");
        REprintf("-------\n");

        REprintf("Time: %g\n", tt);
        REprintf("Node: %i\n", node + 1); /* One based in R */

        REprintf("Current state in node:\n");

        REprintf(" u(length: %i) = {", Nc);
        for (int i = 0; u && i < Nc; i++) {
            REprintf("%i", u[i]);
            if (i < (Nc - 1))
                REprintf(", ");
        }
        REprintf("}\n");

        REprintf(" v(length: %i) = {", Nd);
        for (int i = 0; v && i < Nd; i++) {
            REprintf("%g", v[i]);
            if (i < (Nd - 1))
                REprintf(", ");
        }
        REprintf("}\n");

        REprintf(" ldata(length: %i) = {", Nld);
        for (int i = 0; ldata && i < Nld; i++) {
            REprintf("%g", ldata[i]);
            if (i < (Nld - 1))
                REprintf(", ");
        }
        REprintf("}\n");

        REprintf("Transition: %i\n", transition + 1); /* One based in R */

        if (!R_FINITE(rate) || rate < 0.0)
            REprintf("Rate: %g\n", rate);

        REprintf("\n");

        R_FlushConsole();
    }
}
