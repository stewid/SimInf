/*
 *  siminf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015  Stefan Engblom
 *  Copyright (C) 2015  Stefan Widgren
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

#include <math.h>
#include <string.h>

#include <gsl/gsl_randist.h>

#include "siminf.h"
#include "events.h"

enum {EXIT_EVENT,
      ENTER_EVENT,
      INTERNAL_TRANSFER_EVENT,
      EXTERNAL_TRANSFER_EVENT};

/* Maximum number of individuals to sample from */
#define MAX_INDIVIDUALS 10000
static __thread int kind[MAX_INDIVIDUALS];
static __thread int kind_dest[MAX_INDIVIDUALS];

/**
 * Sample individuals from a node
 *
 * Individuals are sampled from the hidden states determined by select.
 *
 * @param irE Array where irE[k] is the row of E[k]
 * @param jcE jcE[k], index to data of first non-zero element in row k
 * @param Nc Number of compartments in each node.
 * @param state The state vector with number of individuals in each
 *  compartment at each node. The current state in each node is offset
 *  by node * Nc.
 * @param node The node of to sample.
 * @param select Column j in the Select matrix that determines the
 *  hidden states to sample from.
 * @param n The number of individuals to sample. n >= 0.
 * @param proportion If n equals zero, then the number of individuals
 *  to sample is calculated by summing the number of individuals in
 *  the hidden states determined by select and multiplying with the
 *  proportion. 0 <= proportion <= 1.
 * @param individuals The result of the sampling is stored in
 *  the individuals vector.
 * @return 0 on succes or 1 on failure.
 */
static int
sample_select(const size_t *irE,
              const size_t *jcE,
              const size_t Nc,
              const int *state,
              const int node,
              const int select,
              int n,
              double proportion,
              int *inividuals,
              const gsl_rng *rng)
{
    int i, Nstates, Nindividuals = 0, Nkinds = 0;

    /* Clear vector with number of sampled individuals */
    memset(inividuals, 0, Nc * sizeof(int));

    /* 1) Count number of states with individuals */
    /* 2) Count total number of individuals       */
    for (i = jcE[select]; i < jcE[select + 1]; i++) {
        int nk = state[node * Nc + irE[i]];
        if (nk > 0)
            Nkinds++;
        Nindividuals += nk;
    }

    /* Number of hidden states */
    Nstates = jcE[select + 1] - jcE[select];

    /* If n == 0, use the proportion of Nindividuals, else use n as */
    /* the number of individuals to sample                          */
    if (n == 0)
        n = round(proportion * Nindividuals);

    /* Error checking. */
    if (Nstates <= 0        /* No states to sample from, we shouldn't be here. */
        || n > Nindividuals /* Can not sample this number of individuals       */
        || n < 0)           /* Can not sample negative number of individuals.  */
        return 1;

    /* Handle cases that require no random sampling */
    if (n == 0) {
        /* We are done */
        return 0;
    } else if (Nindividuals == n) {
        /* Include all individuals */
        for (i = jcE[select]; i < jcE[select + 1]; i++)
            inividuals[irE[i]] = state[node * Nc + irE[i]];
        return 0;
    } else if (Nstates == 1) {
        /* Only individuals from one state to select from. */
        inividuals[irE[jcE[select]]] = n;
        return 0;
    } else if (Nkinds == 1) {
        /* All individuals to choose from in one state */
        for (i = jcE[select]; i < jcE[select + 1]; i++) {
            if (state[node * Nc + irE[i]] > 0) {
                inividuals[irE[i]] = n;
                break;
            }
        }
        return 0;
    }

    /* Handle cases that require random sampling */
    if (Nstates == 2) {
        /* Sample from the hypergeometric distribution */
        i = jcE[select];
        inividuals[irE[i]] = gsl_ran_hypergeometric(
            rng,
            state[node * Nc + irE[i]],
            state[node * Nc + irE[i+1]],
            n);
        inividuals[irE[i+1]] = n - inividuals[irE[i]];
    } else {
        /* Randomly choose n individuals from a vector of
         * Nindividudals in Nstates */
        int j;

        /* Intialize and populate kind vector */
        if (Nindividuals > MAX_INDIVIDUALS)
            return 1;
        for (i = jcE[select], j = 0; i < jcE[select + 1]; i++) {
            int k, nk, l;

            k  = irE[i];               /* The kind  */
            nk = state[node * Nc + k]; /* N of kind */

            /* Set kind 'k' for 'nk' individuals */
            for (l = 0; l < nk; l++)
                kind[j++] = k;
        }

        /* Randomly choose n individuals from kind vector */
        gsl_ran_choose(rng, kind_dest, n, kind, Nindividuals, sizeof(int));

        /* Count kind of the choosen individuals */
        for (i = 0; i < n; i++)
            inividuals[kind_dest[i]]++;
    }

    return 0;
}

/**
 * Handle exit events
 *
 * Exit events are events that remove individuals from a node.
 *
 * @ingroup events
 * @param irE Array where irE[k] is the row of E[k]
 * @param jcE jcE[k], index to data of first non-zero element in row k
 * @param prE Value of item (i, j) in E.
 * @param Nc Number of compartments in each node.
 * @param Nobs Number of observable states.
 * @param state The state vector with number of individuals in each
 *  compartment at each node. The current state in each node is offset
 *  by node * Nc.
 * @param node The source node of the event.
 * @param dest The dest node of the event.
 * @param select Column j in the Select matrix that determines the
 *  hidden states to sample from.
 * @param n The number of individuals affected by the event. n >= 0.
 * @param proportion If n equals zero, then the number of individuals
 *  affected by the event is calculated by summing the number of
 *  individuals in the hidden states determined by select and
 *  multiplying with the proportion. 0 <= proportion <= 1.
 * @param individuals The result of the sampling is stored in the
 *  individuals vector. Passed as function argument to handle
 *  parallellization.
 * @param rng Random number generator.
 * @return 0 on succes or 1 on failure.
 */
int event_exit(
    const size_t *irE,
    const size_t *jcE,
    const int *prE,
    const size_t Nc,
    const int Nobs,
    int *state,
    const int node,
    const int dest,
    const int select,
    const int n,
    const double proportion,
    int *inividuals,
    const gsl_rng *rng)
{
    int i;

    if (sample_select(irE, jcE, Nc, state, node, select, n, proportion, inividuals, rng))
        return 1;

    for (i = jcE[select]; i < jcE[select + 1]; i++) {
        state[node * Nc + irE[i]] -= inividuals[irE[i]];
        if (state[node * Nc + irE[i]] < 0)
            return SIMINF_ERR_NEGATIVE_STATE;
    }

    return 0;
}

/**
 * Handle enter events
 *
 * Enter events are events that introduce new individuals into a
 * node. All individuals enter the first state for the age category.
 *
 * @ingroup events
 * @param irE Array where irE[k] is the row of E[k]
 * @param jcE jcE[k], index to data of first non-zero element in row k
 * @param prE Value of item (i, j) in E.
 * @param Nc Number of compartments in each node.
 * @param Nobs Number of observable states.
 * @param state The state vector with number of individuals in each
 *  compartment at each node. The current state in each node is offset
 *  by node * Nc.
 * @param node The source node of the event.
 * @param dest The dest node of the event.
 * @param select Column j in the Select matrix that determines the
 *  hidden states to sample from.
 * @param n The number of individuals affected by the event. n >= 0.
 * @param proportion If n equals zero, then the number of individuals
 *  affected by the event is calculated by summing the number of
 *  individuals in the hidden states determined by select and
 *  multiplying with the proportion. 0 <= proportion <= 1.
 * @param individuals The result of the sampling is stored in the
 *  individuals vector. Passed as function argument to handle
 *  parallellization.
 * @param rng Random number generator.
 * @return 0 on succes or 1 on failure.
 */
int event_enter(
    const size_t *irE,
    const size_t *jcE,
    const int *prE,
    const size_t Nc,
    const int Nobs,
    int *state,
    const int node,
    const int dest,
    const int select,
    const int n,
    const double proportion,
    int *inividuals,
    const gsl_rng *rng)
{
    int j = Nobs + select;

    /* :NOTE: All individuals enter first non-zero compartment, i.e. a
     * non-zero entry in element i of select column j */
    if (jcE[j] < jcE[j + 1]) {
        state[node * Nc + irE[jcE[j]]] += n;
        if (state[node * Nc + irE[jcE[j]]] < 0)
            return SIMINF_ERR_NEGATIVE_STATE;
    }

    return 0;
}

/**
 * Handle internal transfer events
 *
 * Internal transfer events are events that change states of
 * individuals whithin one node i.e. ageing n individuals from age_1
 * to age_2.
 *
 * @ingroup events
 * @param irE Array where irE[k] is the row of E[k]
 * @param jcE jcE[k], index to data of first non-zero element in row k
 * @param prE Value of item (i, j) in E.
 * @param Nc Number of compartments in each node.
 * @param Nobs Number of observable states.
 * @param state The state vector with number of individuals in each
 *  compartment at each node. The current state in each node is offset
 *  by node * Nc.
 * @param node The source node of the event.
 * @param dest The dest node of the event.
 * @param select Column j in the Select matrix that determines the
 *  hidden states to sample from.
 * @param n The number of individuals affected by the event. n >= 0.
 * @param proportion If n equals zero, then the number of individuals
 *  affected by the event is calculated by summing the number of
 *  individuals in the hidden states determined by select and
 *  multiplying with the proportion. 0 <= proportion <= 1.
 * @param individuals The result of the sampling is stored in the
 *  individuals vector. Passed as function argument to handle
 *  parallellization.
 * @param rng Random number generator.
 * @return 0 on succes or 1 on failure.
 */
int event_internal_transfer(
    const size_t *irE,
    const size_t *jcE,
    const int *prE,
    const size_t Nc,
    const int Nobs,
    int *state,
    const int node,
    const int dest,
    const int select,
    const int n,
    const double proportion,
    int *inividuals,
    const gsl_rng *rng)
{
    int i, j;

    j = Nobs * INTERNAL_TRANSFER_EVENT + select;
    if (sample_select(irE, jcE, Nc, state, node, j, n, proportion, inividuals, rng))
        return 1;

    for (i = jcE[j]; i < jcE[j + 1]; i++) {
        state[node * Nc + irE[i] + prE[i]] += inividuals[irE[i]];
        if (state[node * Nc + irE[i] + prE[i]] < 0)
            return SIMINF_ERR_NEGATIVE_STATE;
        state[node * Nc + irE[i]] -= inividuals[irE[i]];
        if (state[node * Nc + irE[i]] < 0)
            return SIMINF_ERR_NEGATIVE_STATE;
    }

    return 0;
}

/**
 * Handle external transfer events
 *
 * External transfer events are events that move individuals from one
 * node to another node but keep individuals in the same states
 * i.e. moving n individuals from hidden states of age_1 in node A to
 * the same states of age_1 in node B.
 *
 * @ingroup events
 * @param irE Array where irE[k] is the row of E[k]
 * @param jcE jcE[k], index to data of first non-zero element in row k
 * @param prE Value of item (i, j) in E.
 * @param Nc Number of compartments in each node.
 * @param Nobs Number of observable states.
 * @param state The state vector with number of individuals in each
 *  compartment at each node. The current state in each node is offset
 *  by node * Nc.
 * @param node The source node of the event.
 * @param dest The dest node of the event.
 * @param select Column j in the Select matrix that determines the
 *  hidden states to sample from.
 * @param n The number of individuals affected by the event. n >= 0.
 * @param proportion If n equals zero, then the number of individuals
 *  affected by the event is calculated by summing the number of
 *  individuals in the hidden states determined by select and
 *  multiplying with the proportion. 0 <= proportion <= 1.
 * @param individuals The result of the sampling is stored in the
 *  individuals vector. Passed as function argument to handle
 *  parallellization.
 * @param rng Random number generator.
 * @return 0 on succes or 1 on failure.
 */
int event_external_transfer(
    const size_t *irE,
    const size_t *jcE,
    const int *prE,
    const size_t Nc,
    const int Nobs,
    int *state,
    const int node,
    const int dest,
    const int select,
    const int n,
    const double proportion,
    int *inividuals,
    const gsl_rng *rng)
{
    int i, j;

    j = Nobs * EXTERNAL_TRANSFER_EVENT + select;
    if (sample_select(irE, jcE, Nc, state, node, j, n, proportion, inividuals, rng))
        return 1;

    for (i = jcE[j]; i < jcE[j + 1]; i++) {
        state[dest * Nc + irE[i]] += inividuals[irE[i]];
        if (state[dest * Nc + irE[i]] < 0)
            return SIMINF_ERR_NEGATIVE_STATE;
        state[node * Nc + irE[i]] -= inividuals[irE[i]];
        if (state[node * Nc + irE[i]] < 0)
            return SIMINF_ERR_NEGATIVE_STATE;
    }

    return 0;
}
