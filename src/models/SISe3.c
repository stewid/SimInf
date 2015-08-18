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

#include "SISe3.h"

/* Offset in compartment state vector */
enum {S_1,
      I_1,
      S_2,
      I_2,
      S_3,
      I_3};

/* Offset in model state vector */
enum {PHI};

/* Offsets in node local data (ldata) to parameters in the model */
enum {UPSILON_1,
      UPSILON_2,
      UPSILON_3,
      GAMMA_1,
      GAMMA_2,
      GAMMA_3,
      ALPHA,
      BETA_Q1,
      BETA_Q2,
      BETA_Q3,
      BETA_Q4,
      EPSILON};

/**
 * In age category 1; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param ldata The local data vector for the node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_1_to_I_1(
    const int *u,
    const double *v,
    const double *ldata,
    double t,
    int sd)
{
    return ldata[UPSILON_1] * v[PHI] * u[S_1];
}

/**
 * In age category 2; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param ldata The local data vector for the node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_2_to_I_2(
    const int *u,
    const double *v,
    const double *ldata,
    double t,
    int sd)
{
    return ldata[UPSILON_2] * v[PHI] * u[S_2];
}

/**
 *  In age category 3; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param ldata The local data vector for the node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_3_to_I_3(
    const int *u,
    const double *v,
    const double *ldata,
    double t,
    int sd)
{
    return ldata[UPSILON_3] * v[PHI] * u[S_3];
}

/**
 *  In age category 1; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param ldata The local data vector for the node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_I_1_to_S_1(
    const int *u,
    const double *v,
    const double *ldata,
    double t,
    int sd)
{
    return ldata[GAMMA_1] * u[I_1];
}

/**
 * In age category 2; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param ldata The local data vector for the node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_I_2_to_S_2(
    const int *u,
    const double *v,
    const double *ldata,
    double t,
    int sd)
{
    return ldata[GAMMA_2] * u[I_2];
}

/**
 * In age category 3; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param ldata The local data vector for the node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity
 */
double SISe3_I_3_to_S_3(
    const int *u,
    const double *v,
    const double *ldata,
    double t,
    int sd)
{
    return ldata[GAMMA_3] * u[I_3];
}

/**
 * Update environmental infectious pressure phi
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param ldata The local data vector for the node.
 * @param node The node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return 1 if needs update, else 0.
 */
int SISe3_post_time_step(
    const int *u,
    double *v,
    const double *ldata,
    int node,
    double t,
    int sd)
{
    const int days_in_year = 365;
    const int days_in_quarter = 91;

    double S_n, I_n;
    double tmp = v[PHI];

    S_n = u[S_1] + u[S_2] + u[S_3];
    I_n = u[I_1] + u[I_2] + u[I_3];

    /* Time dependent beta for each quarter of the year. Forward Euler step. */
    switch (((int)t % days_in_year) / days_in_quarter) {
    case 0:
        v[PHI] *= (1.0 - ldata[BETA_Q1]);
        break;
    case 1:
        v[PHI] *= (1.0 - ldata[BETA_Q2]);
        break;
    case 2:
        v[PHI] *= (1.0 - ldata[BETA_Q3]);
        break;
    default:
        v[PHI] *= (1.0 - ldata[BETA_Q4]);
        break;
    }

    if ((I_n + S_n) > 0.0)
        v[PHI] += ldata[ALPHA] * I_n / (I_n + S_n) + ldata[EPSILON];
    else
        v[PHI] += ldata[EPSILON];

    /* 1 if needs update */
    return tmp != v[PHI];
}

/**
 * Run simulation for the SISe3 model
 *
 * This function is called from R with '.Call'
 * @param model The SISe3 model
 * @param threads Number of threads
 * @param seed Random number seed.
 * @return S4 class SISe3 with the simulated trajectory in U
 */
SEXP SISe3_run(SEXP model, SEXP threads, SEXP seed)
{
    int err = 0;
    SEXP result, class_name;
    PropensityFun t_fun[] = {&SISe3_S_1_to_I_1, &SISe3_I_1_to_S_1,
                             &SISe3_S_2_to_I_2, &SISe3_I_2_to_S_2,
                             &SISe3_S_3_to_I_3, &SISe3_I_3_to_S_3};

    if (R_NilValue == model || S4SXP != TYPEOF(model))
        Rf_error("Invalid SISe3 model");

    class_name = getAttrib(model, R_ClassSymbol);
    if (strcmp(CHAR(STRING_ELT(class_name, 0)), "SISe3") != 0)
        Rf_error("Invalid SISe3 model: %s", CHAR(STRING_ELT(class_name, 0)));

    result = PROTECT(duplicate(model));

    err = siminf_run(result, threads, seed, t_fun, &SISe3_post_time_step);

    UNPROTECT(1);

    if (err)
        siminf_error(err);

    return result;
}
