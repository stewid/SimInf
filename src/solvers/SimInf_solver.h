/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2017 - 2018  Robin Eriksson
 *  Copyright (C) 2015 - 2018  Stefan Engblom
 *  Copyright (C) 2015 - 2018  Stefan Widgren
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

#ifndef INCLUDE_SIMINF_SOLVER_H
#define INCLUDE_SIMINF_SOLVER_H

#include <gsl/gsl_rng.h>

#include "kvec.h"
#include "SimInf.h"

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

/* Structure to hold data/arguments to a SimInf solver.
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
*/
typedef struct SimInf_solver_args
{
    /* Initial state vector u0. Integer (Nc X Nn). Gives the initial
     * number of individuals in each compartment in every node. */
    const int *u0;

    /* Initial continuous state vector v0. Double (Nd X Nn). Gives the
     * initial value of the continuous state variables in every
     * node. */
    const double *v0;

    /* Dependency graph. irG[k] is the row of G[k]. */
    const int *irG;

    /* Dependency graph. Index to data of first non-zero element in
     * row k. */
    const int *jcG;

    /* State-change matrix. irS[k] is the row of S[k]. */
    const int *irS;

    /* State-change matrix. Index to data of first non-zero element in
     * row k. */
    const int *jcS;

    /* State-change matrix. Value of item (i, j) in S. */
    const int *prS;

    /* Double vector. Output times. tspan[0] is the start time and
     * tspan[length(tspan)-1] is the stop time. */
    const double *tspan;

    /* Number of sampling points in time. */
    int tlen;

    /* If U is non-NULL, the solution is written to a dense matrix.
     * The output is a matrix U ((Nn * Nc) X length(tspan)). U(:,j)
     * contains the state of the system at tspan(j). */
    int *U;

    /* If U is NULL, the solution is written to a sparse matrix
     * U_sparse. irU[k] is the row of U_sparse[k]. */
    const int *irU;

    /* If U is NULL, the solution is written to a sparse matrix
     * U_sparse. Index to data of first non-zero element in row k. */
    const int *jcU;

    /* If U is NULL, the solution is written to a sparse matrix
     * U_sparse. Value of item (i, j) in U_sparse. */
    double *prU;

    /* If V is non-NULL, the solution is written to a dense matrix.
     * The continuous state output is a matrix V ((Nn * Nd) X
     * length(tspan)).  V(:,j) contains the continuous state of the
     * system at tspan(j). */
    double *V;

    /* If V is NULL, the solution is written to a sparse matrix
     * V_sparse. irV[k] is the row of V_sparse[k]. */
    const int *irV;

    /* If V is NULL, the solution is written to a sparse matrix
     * V_sparse. Index to data of first non-zero element in row k. */
    const int *jcV;

    /* If V is NULL, the solution is written to a sparse matrix
     * V_sparse. Value of item (i, j) in V_sparse. */
    double *prV;

    /* Double matrix (Nld X Nn). Generalized data matrix, data(:,j)
     * gives a local data vector for node #j. */
    const double *ldata;

    /* The global data vector. */
    const double *gdata;

    /* Number of nodes. */
    int Nn;

    /* Number of compartments in each node. */
    int Nc;

    /* Number of different transitions. */
    int Nt;

    /* Number of continuous state variables. */
    int Nd;

    /* Length of the local data vector 'ldata' for each node. The
     * 'ldata' vector is sent to propensities and the post time step
     * function. */
    int Nld;

    /* Select matrix for events. irE[k] is the row of E[k]. */
    const int *irE;

    /* Select matrix for events. Index to data of first non-zero
     * element in row k. */
    const int *jcE;

    /* Shift matrix for internal and external transfer events. */
    const int *N;

    /* Number of events. */
    int len;

    /* The type of event i. */
    const int *event;

    /* The time of event i. */
    const int *time;

    /* The source node of event i. */
    const int *node;

    /* The dest node of event i. */
    const int *dest;

    /* The number of individuals in the scheduled event. n[i] >= 0. */
    const int *n;

    /* If n[i] equals zero, then the number of individuals to sample
     * is calculated by summing the number of individuals in the
     * states determined by select[i] and multiplying with the
     * proportion. 0 <= p[i] <= 1. */
    const double *proportion;

    /* Column j in the event matrix E that determines the states to
     * sample from. */
    const int *select;

    /* Column j in the shift matrix S that determines the shift of the
     * internal and external transfer event. */
    const int *shift;

    /* Number of threads to use during simulation. */
    int Nthread;

    /* Random number seed. */
    unsigned long int seed;

    /* Vector of function pointers to transition rate functions. */
    TRFun *tr_fun;

    /* Function pointer to callback after each time step e.g. to
     * update the infectious pressure. */
    PTSFun pts_fun;
} SimInf_solver_args;

/**
 * Structure with data for a scheduled event.
 */
typedef struct SimInf_scheduled_event
{
    int event;         /**< The type of the event. */
    int time;          /**< The time for the event. */
    int node;          /**< The source node of the event. */
    int dest;          /**< The dest node of the event. */
    int n;             /**< The number of individuals in the scheduled
                        *   event. n >= 0. */
    double proportion; /**< If n equals zero, then the number of
                        *   individuals to sample is calculated by
                        *   summing the number of individuals in the
                        *   states determined by select and
                        *   multiplying with the proportion.  0 <=
                        *   proportion <= 1. */
    int select;        /**< Column j in the event matrix that
                        *   determines the states to sample from. */
    int shift;         /**< Column j in the shift matrix that
                        *   determines the shift of the internal
                        *   and external transfer event. */
} SimInf_scheduled_event;

typedef kvec_t(SimInf_scheduled_event) SimInf_events_t;

/**
 * Structure with data to process scheduled events.
 */
typedef struct SimInf_scheduled_events
{
    /*** Constants ***/
    int Nthread;          /**< Number of threads. */

    /*** Matrices to process events ***/
    const int *irE;       /**< Select matrix for events. irE[k] is the
                           *   row of E[k]. */
    const int *jcE;       /**< Select matrix for events. Index to data
                           *   of first non-zero element in row k. */
    const int *N;         /**< Shift matrix for internal and external
                           *   transfer events. */

    /*** Scheduled events ***/
    SimInf_events_t events; /**< Events to process. */
    int events_index;       /**< Index to the next event to
                             *   process. */

    /*** Vectors for sampling individuals ***/
    int *individuals;     /**< Vector to store the result of the
                           *   sampling during scheduled events
                           *   processing. */
    int *u_tmp;           /**< Temporary vector with the compartment
                           *   state in a node when sampling
                           *   individuals for scheduled events. */
    gsl_rng *rng;         /**< The random number generator for
                           *   sampling. */
} SimInf_scheduled_events;

/**
 * Structure to hold thread specific data/arguments for simulation.
 */
typedef struct SimInf_compartment_model
{
    /*** Constants ***/
    int Nthread; /**< Number of threads. */
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

    /*** Callbacks ***/
    TRFun *tr_fun;  /**< Vector of function pointers to
                     *   transition rate functions */
    PTSFun pts_fun; /**< Callback after each time step */

    /*** Keep track of time ***/
    double tt;           /**< The global time. */
    double next_unit_of_time; /**< The global time of next unit of
                               * time. */
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
    int *update_node; /**< Vector of length Nn used to indicate nodes
                       *   for update. */

    double *sum_t_rate; /**< Vector of length Nn with the sum of
                         *   propensities in every node. */
    double *t_rate;     /**< Transition rate matrix (Nt X Nn) with all
                         *   propensities for state transitions. */
    double *t_time;     /**< Time for next event (transition) in each
                         *   node. */
    int error;          /**< The error state of the thread. 0 if
                         *   ok. */
} SimInf_compartment_model;

int SimInf_compartment_model_create(
    SimInf_compartment_model **out, SimInf_solver_args *args);

void SimInf_compartment_model_free(
    SimInf_compartment_model *model);

int SimInf_scheduled_events_create(
    SimInf_scheduled_events **out, SimInf_solver_args *args, gsl_rng *rng);

void SimInf_scheduled_events_free(
    SimInf_scheduled_events *events);

void SimInf_process_events(
    SimInf_compartment_model *model,
    SimInf_scheduled_events *events,
    int process_E2);

void SimInf_store_solution_sparse(SimInf_compartment_model *model);

#endif
