/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 -- 2026 Stefan Widgren
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

#include "SimInf_internal.h"
#include <R_ext/Visibility.h>

typedef kvec_t(int) kvec_t_int;

typedef struct SimInf_raster_model
{
    /*** Data vectors for propensities ***/
    double *sum_rate;      /**< Vector of length Ncells with the sum
                            *   of cell and node propensities in every
                            *   cell. */
    double *sum_node_rate; /**< Vector of length Nnodes with the sum
                            *   of propensities in every node. */
    double *node_rate;     /**< Transition rate matrix (Nt X Nnodes)
                            *   with all propensities for state
                            *   transitions in nodes. */
    double *sum_cell_rate; /**< Vector of length Ncells with the sum
                            *   of propensities in every cell. */
    double *cell_rate;     /**< Transition rate matrix (Nt X Ncells)
                            *   with all propensities for state
                            *   transitions in cells. */

    /*** Binary (min)heap. ***/
    int *cells;
    int *heap;

    /*** Keep track of nodes ***/
    kvec_t_int *nodes;

    /*** Callbacks ***/
    TRRasterFun *tr_fun;  /**< Vector of function pointers to
                           *   transition rate functions. */
    PTSFun pts_fun;       /**< Callback after each time unit step. */

    /*** Keep track of time ***/
    const double *tspan; /**< Output times. tspan[0] is the start time
                          *   and tspan[length(tspan)-1] is the stop
                          *   time.*/
    int tlen;            /**< Number of sampling points in time. */
    int t_it;            /**< Index to next time in tspan */
    double *cell_time;   /**< Time for next event in each cell. */

    /*** Dependency graph ***/
    const int *irG; /**< Dependency graph. irG[k] is the row of
                     *   G[k]. */
    const int *jcG; /**< Dependency graph. Index to data of first
                     *   non-zero element in row k. */

    /*** Data vectors for the nodes ***/
    int *u;           /**< Vector with the number of individuals in
                       *   each compartment in each node in the
                       *   thread. */
    int *U;           /**< If the solution is written to a dense
                       *   matrix the compartment output is a matrix U
                       *   ((Nnodes * Nc) X length(tspan)). U(:,j)
                       *   contains the state of the system at
                       *   tspan(j). */
    const int *irU;   /**< If the solution is written to a sparse
                       *   matrix, irU[k] is the row of U[k]. */
    const int *jcU;   /**< If the solution is written to a sparse
                       *   matrix, index to data of first non-zero
                       *   element in row k. */
    double    *prU;   /**< If the solution is written to a sparse
                       *   matrix, value of item (i, j) in U. */
    const int *irS;   /**< Node state-change matrix. irS[k] is the row
                       *   of S[k]. */
    const int *jcS;   /**< Node state-change matrix. Index to data of
                       *   first non-zero element in row k. */
    const int *prS;   /**< Node state-change matrix. Value of item (i,
                       *   j) in S. */

    const double *ldata; /**< Matrix (Nldata X Nnodes). ldata(:,j)
                          *   gives a local data vector for node
                          *   #j. */
    const double *gdata; /**< The global data vector. */

    /*** Data vectors for the cells ***/
    int *cell_u;  /**< Vector with the count of each compartment in
                   *   each cell. */

    /*** Sparse matrices ***/
    const int *cell_irS; /**< Cell state-change matrix. cell_irS[k] is
                          *   the row of cell_S[k]. */
    const int *cell_jcS; /**< Cell state-change matrix. Index to data
                          *   of first non-zero element in row k. */
    const int *cell_prS; /**< Cell state-change matrix. Value of item
                          *   (i, j) in cell_S. */

    /*** Transition type ***/
    const int *tr_type;  /**< Keep track of if a transition happens on
                          *   a cell, in a node or if it is a
                          *   movement. */

    gsl_rng *rng;        /**< The random number generator. */

    /*** Constants ***/
    int Nnodes;     /**< Number of nodes. */
    int Ncells;     /**< Number of cells. */
    int Nc;         /**< Number of compartments in each node. */
    int cell_Nc;    /**< Number of compartments in each cell. */
    int cell_i;     /**< Index to the cell compartment in each
                     *   node. */
    int Nt;         /**< Total number of different transitions. */
    int Nldata;     /**< Length of the local data vector 'ldata' for
                     *   each node. The 'ldata' vector is sent to the
                     *   propensity functions. */
} SimInf_raster_model;

/**
 * Initialize and run the SimInf raster solver.
 *
 * @param args Structure with data for the solver.
 * @return 0 if Ok, else error code.
 */
attribute_hidden int
SimInf_run_solver_raster(
    SimInf_solver_args *args)
{
    /* Not yet implemented. */
    return 0;
}
