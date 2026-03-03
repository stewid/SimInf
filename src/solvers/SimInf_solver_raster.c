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
    double *prV;      /**< If the solution is written to a sparse
                       *   matrix, value of item (i, j) in V. */
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
    int nrow;       /**< Number of rows. */
    int ncol;       /**< Number of cols. */
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
 * Free allocated memory for a raster model.
 *
 * @param model the data structure to free.
 */
static void
SimInf_raster_model_free(
    SimInf_raster_model *model)
{
    if (model) {
        /* Free data vectors for propensities. */
        free(model->sum_rate);
        model->sum_rate = NULL;
        free(model->sum_node_rate);
        model->sum_node_rate = NULL;
        free(model->node_rate);
        model->node_rate = NULL;
        free(model->sum_cell_rate);
        model->sum_cell_rate = NULL;
        free(model->cell_rate);
        model->cell_rate = NULL;

        /* Free data for binary (min)heap. */
        free(model->cells);
        model->cells = NULL;
        free(model->heap);
        model->heap = NULL;

        /* Free data to keep track of nodes. */
        if (model->nodes) {
            const int ncells = model->nrow * model->ncol;
            for (int i = 0; i < ncells; i++)
                kv_destroy(model->nodes[i]);
            free(model->nodes);
            model->nodes = NULL;
        }

        /* Free data to keep track of time. */
        free(model->cell_time);
        model->cell_time = NULL;

        /* Free data vectors to keep track of the states in the
         * nodes. */
        free(model->u);
        model->u = NULL;
        free(model->cell_u);
        model->cell_u = NULL;
        free(model->v);
        model->v = NULL;
        free(model->v_new);
        model->v_new = NULL;

        /* Free data for the random number generator. */
        gsl_rng_free(model->rng);
        model->rng = NULL;

        free(model);
    }
}

/**
 * Create and initialize data to process a raster model. The generated
 * data structure must be freed by the user.
 *
 * @param out the resulting data structure.
 * @param args structure with data for the solver.
 * @param rng random number generator
 * @return 0 or an error code
 */
static int
SimInf_raster_model_create(
    SimInf_raster_model **out,
    SimInf_solver_args const *args,
    gsl_rng *rng)
{
    int err = SIMINF_ERR_ALLOC_MEMORY_BUFFER;
    SimInf_raster_model *model = NULL;

    model = calloc(1, sizeof(SimInf_raster_model));
    if (!model)
        goto on_error;

    /* Callbacks */
    model->tr_fun = args->tr_raster_fun;
    model->pts_fun = args->pts_fun;

    /* Constants. */
    model->Nnodes = args->Nn;
    model->Nc = args->Nc;
    model->nrow = args->nrow;
    model->ncol = args->ncol;
    if (model->nrow < 1 || model->ncol < 1)
        goto on_error;
    model->Nt = args->Nt;
    model->Nldata = args->Nld;

    /* Keep track of time. */
    model->tspan = args->tspan;
    model->tlen = args->tlen;

    /* Number of compartments in each cell. */
    model->cell_Nc = args->cell_Nc;
    if (model->cell_Nc < 0)
        goto on_error;

    /* Index to the cell compartment in each node. */
    model->cell_i = args->cell_i;
    if (model->cell_i < 0 || model->cell_i >= args->Nc)
        goto on_error;

    /* Nodes */
    const ptrdiff_t n_cells = (ptrdiff_t) model->nrow * (ptrdiff_t) model->ncol;
    model->nodes = calloc(n_cells, sizeof(kvec_t_int));
    if (!model->nodes)
        goto on_error;

    /* Dependency graph. */
    model->irG = args->irG;
    model->jcG = args->jcG;

    /* Binary (min)heap. */
    model->cells = malloc(n_cells * sizeof(int));
    if (!model->cells)
        goto on_error;

    model->heap = malloc(n_cells * sizeof(int));
    if (!model->heap)
        goto on_error;

    model->cell_time = malloc(n_cells * sizeof(double));
    if (!model->cell_time)
        goto on_error;

    /* Data vectors for the nodes. */
    if (args->U) {
        model->U = args->U;
    } else {
        model->irU = args->irU;
        model->jcU = args->jcU;
        model->prU = args->prU;
    }
    if (args->V) {
        model->V = args->V;
    } else {
        model->irV = args->irV;
        model->jcV = args->jcV;
        model->prV = args->prV;
    }
    model->irS = args->irS;
    model->jcS = args->jcS;
    model->prS = args->prS;

    /* Check that the state-change-matrix doesn't change the value in
     * the cell compartment. */
    for (int i = 0; i < model->Nt; i++) {
        for (int j = model->jcS[i]; j < model->jcS[i + 1]; j++) {
            if (model->irS[j] == model->cell_i) {
                err = SIMINF_ERR_NON_ZERO_CELL_IN_S;
                goto on_error;
            }
        }
    }

    /* State-change matrix for the cell. */
    model->cell_irS = args->cell_irS;
    model->cell_jcS = args->cell_jcS;
    model->cell_prS = args->cell_prS;

    /* Local and global data. */
    model->ldata = args->ldata;
    model->gdata = args->gdata;

    /* Allocate memory for the cell compartment state and set it to
     * the initial state. */
    const ptrdiff_t cell_Nc = args->cell_Nc;
    model->cell_u = calloc(n_cells * cell_Nc, sizeof(int));
    if (!model->cell_u)
        goto on_error;

    /* Allocate memory for the node compartment state and set it to
     * the initial state. */
    const ptrdiff_t Nn = args->Nn;
    const ptrdiff_t Nc = args->Nc;
    model->u = malloc(Nn * Nc * sizeof(int));
    if (!model->u)
        goto on_error;
    memcpy(model->u, args->u0, Nn * Nc * sizeof(int));

    /* Allocate memory to keep track of the continuous state in each
     * node. */
    const ptrdiff_t Nd = args->Nd;
    const ptrdiff_t v_len = Nn * Nd * sizeof(double);
    model->v = malloc(v_len);
    if (!model->v)
        goto on_error;          /* #nocov */
    model->v_new = malloc(v_len);
    if (!model->v_new)
        goto on_error;          /* #nocov */

    /* Set continuous state to the initial state in each node. */
    if (v_len > 0) {
        memcpy(model[0].v, args->v0, v_len);
        memcpy(model[0].v_new, args->v0, v_len);
    }

    /* Create transition rate matrix (Nt X Nnodes) and total rate
     * vector. In node_rate we store all propensities for state
     * transitions, and in sum_node_rate the sum of propensities in
     * every node. */
    model->node_rate = calloc(model->Nt * model->Nnodes, sizeof(double));
    if (!model->node_rate)
        goto on_error;
    model->sum_node_rate = calloc(model->Nnodes, sizeof(double));
    if (!model->sum_node_rate)
        goto on_error;

    /* Create transition rate matrix (Nt X Ncells) and total rate
     * vector. In cell_rate we store all propensities for state
     * transitions, and in sum_cell_rate the sum of propensities in
     * every cell. */
    model->cell_rate = calloc(model->Nt * n_cells, sizeof(double));
    if (!model->cell_rate)
        goto on_error;
    model->sum_cell_rate = calloc(n_cells, sizeof(double));
    if (!model->sum_cell_rate)
        goto on_error;

    model->sum_rate = calloc(n_cells, sizeof(double));
    if (!model->sum_rate)
        goto on_error;

    /* Random number generator */
    model->rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!model->rng)
        goto on_error;
    gsl_rng_set(model->rng, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));

    *out = model;
    return 0;

on_error:
    SimInf_raster_model_free(model);
    return err;
}

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
