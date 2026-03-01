/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
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

#pragma once

#include <R.h>
#include <Rinternals.h>

#ifndef R_PRIdXLEN_T
#  ifdef LONG_VECTOR_SUPPORT
#    define R_PRIdXLEN_T "td"
#  else
#    define R_PRIdXLEN_T "d"
#  endif
#endif

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

#define SIMINF_UNUSED(x) ((void)(x))
#define SIMINF_STR(name) #name
#define SIMINF_CALLDEF(name, n) {SIMINF_STR(name), (DL_FUNC) &name, n}

/* Error constants */
typedef enum {
    SIMINF_ERR_NEGATIVE_STATE       = -1,
    SIMINF_ERR_ALLOC_MEMORY_BUFFER  = -2,
    SIMINF_UNDEFINED_EVENT          = -3,
    SIMINF_INVALID_EDGE_PROBABILITY = -4,
    SIMINF_ERR_V_IS_NOT_FINITE      = -7,
    SIMINF_ERR_SAMPLE_SELECT        = -8,
    SIMINF_ERR_INVALID_MODEL        = -9,
    SIMINF_ERR_V_IS_NEGATIVE        = -10,
    SIMINF_ERR_INVALID_RATE         = -11,
    SIMINF_ERR_UNKNOWN_SOLVER       = -12,
    SIMINF_ERR_DEST_OUT_OF_BOUNDS   = -13,
    SIMINF_ERR_NODE_OUT_OF_BOUNDS   = -14,
    SIMINF_ERR_EVENTS_N             = -15,
    SIMINF_ERR_EVENT_SHIFT          = -16,
    SIMINF_ERR_SHIFT_OUT_OF_BOUNDS  = -17,
    SIMINF_ERR_INVALID_PROPORTION   = -18,
    SIMINF_ERR_AEM_REPLICATED_MODEL = -19
} SimInf_error_code;

/* Forward declaration of the transition rate function for a
 * SimInf_model. */
typedef double (*TRFun)(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t);

/* Forward declaration of the transition rate function for a
 * SimInf_raster_model. */
typedef double (*TRRasterFun)(
    const int *cell,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t);

/* Forward declaration of the post time step callback function. */
typedef int (*PTSFun)(
    double *v_new,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    int node,
    double t);

/* Forward declaration of the function to initiate and run the
 * simulation for a SimInf_model object. */
SEXP SimInf_run(
    SEXP model,
    SEXP solver,
    TRFun *tr_fun,
    PTSFun pts_fun);

/* Forward declaration of the function to initiate and run the
 * simulation for a SimInf_raster_model object. */
SEXP SimInf_raster_run(
    SEXP model,
    TRRasterFun *tr_fun,
    PTSFun pts_fun);

/**
 * Decay of environmental infectious pressure with a forward Euler
 * step.
 *
 * The time dependent beta is divided into four intervals of the year
 * where 0 <= day < 365
 *
 * Case 1: END_1 < END_2 < END_3 < END_4
 * INTERVAL_1 INTERVAL_2     INTERVAL_3     INTERVAL_4     INTERVAL_1
 * [0, END_1) [END_1, END_2) [END_2, END_3) [END_3, END_4) [END_4, 365)
 *
 * Case 2: END_3 < END_4 < END_1 < END_2
 * INTERVAL_3 INTERVAL_4     INTERVAL_1     INTERVAL_2     INTERVAL_3
 * [0, END_3) [END_3, END_4) [END_4, END_1) [END_1, END_2) [END_2, 365)
 *
 * Case 3: END_4 < END_1 < END_2 < END_3
 * INTERVAL_4 INTERVAL_1     INTERVAL_2     INTERVAL_3     INTERVAL_4
 * [0, END_4) [END_4, END_1) [END_1, END_2) [END_2, END_3) [END_3, 365)
 *
 * @param phi The currrent value of the environmental infectious pressure.
 * @param day The day of the year 0 <= day < 365.
 * @param end_t1 The non-inclusive day that ends interval 1.
 * @param end_t2 The non-inclusive day that ends interval 2.
 * @param end_t3 The non-inclusive day that ends interval 3.
 * @param end_t4 The non-inclusive day that ends interval 4.
 * @param beta_t1 The value for beta in interval 1.
 * @param beta_t2 The value for beta in interval 2.
 * @param beta_t3 The value for beta in interval 3.
 * @param beta_t4 The value for beta in interval 4.
 * @return phi * (1.0 - beta) (where beta is the value for the interval)
 */
double SimInf_forward_euler_linear_decay(
    double phi, int day,
    int end_t1, int end_t2, int end_t3, int end_t4,
    double beta_t1, double beta_t2, double beta_t3, double beta_t4);

/**
 * Local spread of the environmental infectious pressure phi among
 * proximal nodes.
 *
 * @param neighbors Spatial coupling between nodes where 'neighbors'
 * is a vector of pairs (index, distance) to neighbor nodes. The pair
 * vector is terminated with an index equal to -1.
 * @param phi Vector with phi in each node
 * @param u The compartment state vector in each node.
 * @param N_i The number of individuals in node i.
 * @param phi_i The environmental infectious pressure phi in node i.
 * @param Nc The number of compartments in each node.
 * @param D The spatial coupling of the environmental infectious
 * pressure phi among proximal nodes.
 * @return The contribution from neighbors to phi in node i
 */
double SimInf_local_spread(
    const double *neighbors, const double *phi,
    const int *u, const double N_i,
    const double phi_i, const int Nc, const double D);
