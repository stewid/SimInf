/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015 - 2016  Stefan Engblom
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

#include "SISe3_sp.h"
#include "siminf_forward_euler_linear_decay.h"

/* Offset in integer compartment state vector */
enum {S_1, I_1, S_2, I_2, S_3, I_3};

/* Offset in real-valued continuous state vector */
enum {PHI};

/* Offsets in node local data (ldata) to parameters in the model */
enum {END_T1, END_T2, END_T3, END_T4, NEIGHBOR};

/* Offsets in global data (gdata) to parameters in the model */
enum {UPSILON_1, UPSILON_2, UPSILON_3, GAMMA_1, GAMMA_2, GAMMA_3,
      ALPHA, BETA_T1, BETA_T2, BETA_T3, BETA_T4, EPSILON, COUPLING};

/**
 * In age category 1; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_sp_S_1_to_I_1(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[UPSILON_1] * v[PHI] * u[S_1];
}

/**
 * In age category 2; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_sp_S_2_to_I_2(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[UPSILON_2] * v[PHI] * u[S_2];
}

/**
 *  In age category 3; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_sp_S_3_to_I_3(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[UPSILON_3] * v[PHI] * u[S_3];
}

/**
 *  In age category 1; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_sp_I_1_to_S_1(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[GAMMA_1] * u[I_1];
}

/**
 * In age category 2; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_sp_I_2_to_S_2(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[GAMMA_2] * u[I_2];
}

/**
 * In age category 3; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity
 */
double SISe3_sp_I_3_to_S_3(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[GAMMA_3] * u[I_3];
}

/**
 * Update environmental infectious pressure phi
 *
 * Decay environmental infectious pressure phi, add contribution from
 * infected individuals and proximity coupling.
 * @param v_new The continuous state vector in the node after the post
 * time step
 * @param u The compartment state vector in the node.
 * @param v The current continuous state vector in the node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param node The node.
 * @param t Current time.
 * @param sd The sub-domain of the node.
 * @return error code (<0), or 1 if node needs to update the
 * transition rates, or 0 when it doesn't need to update the
 * transition rates.
 */
int SISe3_sp_post_time_step(
    double *v_new,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    int node,
    double t,
    int sd)
{
    int i, j;
    const int day = (int)t % 365;
    const double I_n = u[I_1] + u[I_2] + u[I_3];
    const double n = u[S_1] + u[S_2] + u[S_3] + I_n;
    const double phi = v[PHI];
    const double coupling = gdata[COUPLING];

    /* Deterimine the pointer to the continuous state vector in the
     * first node. Use this to find phi at neighbours to the current
     * node. */
    const double *phi_0 = &v[-node];

    /* Time dependent beta in each of the four intervals of the
     * year. Forward Euler step. */
    v_new[PHI] = siminf_forward_euler_linear_decay(
        phi, day,
        ldata[END_T1], ldata[END_T2], ldata[END_T3], ldata[END_T4],
        gdata[BETA_T1], gdata[BETA_T2], gdata[BETA_T3], gdata[BETA_T4]);

    if (n > 0.0)
        v_new[PHI] += gdata[ALPHA] * I_n / n + gdata[EPSILON];
    else
        v_new[PHI] += gdata[EPSILON];

    /* Coupling between neighboring farms. */
    /* i is the offset in local data to the first neighbor. */
    /* j is the neighbor node or -1 to stop.  */
    i = NEIGHBOR;
    j = (int)ldata[i];
    while (j >= 0) {
        v_new[PHI] += (phi_0[j] - phi) * coupling * ldata[i + 1];

        /* Move to next neighbor pair (index, value) */
        i += 2;
        j = (int)ldata[i];
    }

    if (isfinite(v_new[PHI]))
        return phi != v_new[PHI]; /* 1 if needs update */
    return SIMINF_ERR_V_IS_NOT_FINITE;
}

/**
 * Run simulation for the SISe3_sp model
 *
 * This function is called from R with '.Call'
 * @param model The SISe3_sp model
 * @param threads Number of threads
 * @param seed Random number seed.
 * @return S4 class SISe3_sp with the simulated trajectory in U
 */
SEXP SISe3_sp_run(SEXP model, SEXP threads, SEXP seed)
{
    int err = 0;
    SEXP result, class_name;
    PropensityFun t_fun[] = {&SISe3_sp_S_1_to_I_1, &SISe3_sp_I_1_to_S_1,
                             &SISe3_sp_S_2_to_I_2, &SISe3_sp_I_2_to_S_2,
                             &SISe3_sp_S_3_to_I_3, &SISe3_sp_I_3_to_S_3};

    if (R_NilValue == model || S4SXP != TYPEOF(model))
        Rf_error("Invalid SISe3_sp model");

    class_name = getAttrib(model, R_ClassSymbol);
    if (strcmp(CHAR(STRING_ELT(class_name, 0)), "SISe3_sp") != 0)
        Rf_error("Invalid SISe3_sp model: %s", CHAR(STRING_ELT(class_name, 0)));

    result = PROTECT(duplicate(model));

    err = siminf_run(result, threads, seed, t_fun, &SISe3_sp_post_time_step);

    UNPROTECT(1);

    if (err)
        siminf_error(err);

    return result;
}
