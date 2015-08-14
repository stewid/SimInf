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

/* Offset in compartment state vector */
enum {S_1,
      I_1,
      S_2,
      I_2,
      S_3,
      I_3};

/* Offset in model state vector */
enum {PHI};

/* Offsets in data to parameters in the model */
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
 * @param data The data vector for node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_1_to_I_1(
    const int *u,
    const double *v,
    const double *data,
    double t,
    int sd)
{
    return data[UPSILON_1] * v[PHI] * u[S_1];
}

/**
 * In age category 2; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param data The data vector for node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_2_to_I_2(
    const int *u,
    const double *v,
    const double *data,
    double t,
    int sd)
{
    return data[UPSILON_2] * v[PHI] * u[S_2];
}

/**
 *  In age category 3; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param data The data vector for node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_3_to_I_3(
    const int *u,
    const double *v,
    const double *data,
    double t,
    int sd)
{
    return data[UPSILON_3] * v[PHI] * u[S_3];
}

/**
 *  In age category 1; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param data The data vector for node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_I_1_to_S_1(
    const int *u,
    const double *v,
    const double *data,
    double t,
    int sd)
{
    return data[GAMMA_1] * u[I_1];
}

/**
 * In age category 2; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param data The data vector for node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_I_2_to_S_2(
    const int *u,
    const double *v,
    const double *data,
    double t,
    int sd)
{
    return data[GAMMA_2] * u[I_2];
}

/**
 * In age category 3; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param data The data vector for node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity
 */
double SISe3_I_3_to_S_3(
    const int *u,
    const double *v,
    const double *data,
    double t,
    int sd)
{
    return data[GAMMA_3] * u[I_3];
}

/**
 * Update infectious pressure
 *
 * @param u The compartment state vector in node.
 * @param v The model state vector in node.
 * @param data The data vector for node.
 * @param node The node.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return 1 if needs update, else 0.
 */
int SISe3_post_time_step(
    const int *u,
    double *v,
    const double *data,
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
        v[PHI] *= (1.0 - data[BETA_Q1]);
        break;
    case 1:
        v[PHI] *= (1.0 - data[BETA_Q2]);
        break;
    case 2:
        v[PHI] *= (1.0 - data[BETA_Q3]);
        break;
    default:
        v[PHI] *= (1.0 - data[BETA_Q4]);
        break;
    }

    if ((I_n + S_n) > 0.0)
        v[PHI] += data[ALPHA] * I_n / (I_n + S_n) + data[EPSILON];
    else
        v[PHI] += data[EPSILON];

    /* 1 if needs update */
    return tmp != v[PHI];
}
