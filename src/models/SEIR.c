/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
 * Copyright (C) 2015 -- 2020 Stefan Widgren
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

#include "SimInf.h"

/* Offset in integer compartment state vector */
enum {S, E, I, R};

/* Offsets in global data (gdata) to parameters in the model */
enum {BETA, EPSILON, GAMMA};

/**
 * susceptible to exposed: S -> E
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @return propensity.
 */
double SEIR_S_to_E(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t)
{
    const double S_n = u[S];
    const double I_n = u[I];
    const double n = S_n + u[E] + I_n + u[R];

    if (n > 0.0)
        return gdata[BETA] * S_n * I_n / n;
    return 0.0;
}

/**
 * exposed to infected: E -> I
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @return propensity.
 */
double SEIR_E_to_I(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t)
{
    return gdata[EPSILON] * u[E];
}

/**
 *  infected to recovered: I -> R
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @return propensity.
 */
double SEIR_I_to_R(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t)
{
    return gdata[GAMMA] * u[I];
}

/**
 * SEIR post time step
 *
 * @param v_new The continuous state vector in the node after the post
 * time step
 * @param u The compartment state vector in the node.
 * @param v The current continuous state vector in the node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param node The node.
 * @param t The current time.
 * @return error code (<0), or 1 if node needs to update the
 * transition rates, or 0 when it doesn't need to update the
 * transition rates.
 */
int SEIR_post_time_step(
    double *v_new,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    int node,
    double t)
{
    return 0;
}

/**
 * Run simulation with the SEIR model
 *
 * @param model The SIR model.
 * @param threads Number of threads.
 * @param solver The numerical solver.
 * @return The simulated trajectory.
 */
SEXP SEIR_run(SEXP model, SEXP threads, SEXP solver)
{
    TRFun tr_fun[] = {&SEIR_S_to_E, &SEIR_E_to_I, &SEIR_I_to_R};

    return SimInf_run(model, threads, solver, tr_fun, &SEIR_post_time_step);
}
