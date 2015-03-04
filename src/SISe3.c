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

/* Compartments */
enum {S_1,
      I_1,
      S_2,
      I_2,
      S_3,
      I_3};

/* Offsets in data to parameters in the model */
enum {PHI,
      UPSILON_1,
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
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_1_to_I_1(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return data[UPSILON_1] * data[PHI] * x[S_1];
}

/**
 * In age category 2; susceptible to infected: S -> I
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_2_to_I_2(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return data[UPSILON_2] * data[PHI] * x[S_2];
}

/**
 *  In age category 3; susceptible to infected: S -> I
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_3_to_I_3(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return data[UPSILON_3] * data[PHI] * x[S_3];
}

/**
 *  In age category 1; infected to susceptible: I -> S
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_I_1_to_S_1(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return data[GAMMA_1] * x[I_1];
}

/**
 * In age category 2; infected to susceptible: I -> S
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_I_2_to_S_2(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return data[GAMMA_2] * x[I_2];
}

/**
 * In age category 3; infected to susceptible: I -> S
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity
 */
double SISe3_I_3_to_S_3(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return data[GAMMA_3] * x[I_3];
}

/**
 * Update infectious pressure
 *
 * @param x The state vector in node.
 * @param node The node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return 1 if needs update, else 0.
 */
int SISe3_post_time_step(
    const int *x,
    int node,
    double t,
    double *data,
    int sd)
{
    const int days_in_year = 365;
    const int days_in_quarter = 91;

    double S_n, I_n;
    double tmp = data[PHI];

    S_n = x[S_1] + x[S_2] + x[S_3];
    I_n = x[I_1] + x[I_2] + x[I_3];

    /* Time dependent beta for each quarter of the year. Forward Euler step. */
    switch (((int)t % days_in_year) / days_in_quarter) {
    case 0:
        data[PHI] *= (1.0 - data[BETA_Q1]);
        break;
    case 1:
        data[PHI] *= (1.0 - data[BETA_Q2]);
        break;
    case 2:
        data[PHI] *= (1.0 - data[BETA_Q3]);
        break;
    default:
        data[PHI] *= (1.0 - data[BETA_Q4]);
        break;
    }
    data[PHI] += data[ALPHA] * I_n / (I_n + S_n) + data[EPSILON];

    /* 1 if needs update */
    return tmp != data[PHI];
}
