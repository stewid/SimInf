/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2017  Robin Eriksson
 *  Copyright (C) 2015 - 2017  Stefan Engblom
 *  Copyright (C) 2015 - 2017  Stefan Widgren
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

#include "SimInf.h"

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

#endif
