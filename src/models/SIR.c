/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015 - 2016  Stefan Widgren
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

#include "siminf.h"

/* Offset in integer compartment state vector */
enum {S, I, R};

/* Offsets in global data (gdata) to parameters in the model */
enum {BETA, GAMMA};

/**
 * susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @return propensity.
 */
double SIR_S_to_I(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t)
{
    const double S_n = u[S];
    const double I_n = u[I];

    return (gdata[BETA] * S_n * I_n) / (S_n + I_n + u[R]);
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
double SIR_I_to_R(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t)
{
    return gdata[GAMMA] * u[I];
}

/**
 * SIR post time step
 *
 * @param v_new The continuous state vector in the node after the post
 * time step
 * @param u The compartment state vector in the node.
 * @param v The current continuous state vector in the node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param node The node.
 * @param t The current time.
 * @param rng The random number generator.
 * @return error code (<0), or 1 if node needs to update the
 * transition rates, or 0 when it doesn't need to update the
 * transition rates.
 */
int SIR_post_time_step(
    double *v_new,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    int node,
    double t,
    gsl_rng *rng)
{
    return 0;
}

/**
 * Run simulation with the SIR model
 *
 * @param model The SIR model.
 * @param threads Number of threads.
 * @param seed Random number seed.
 * @return The simulated trajectory.
 */
SEXP SIR_run(SEXP model, SEXP threads, SEXP seed)
{
    TRFun tr_fun[] = {&SIR_S_to_I, &SIR_I_to_R};

    return siminf_run(model, threads, seed, tr_fun, &SIR_post_time_step);
}
