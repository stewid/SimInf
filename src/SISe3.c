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
enum {S_age_1,
      I_age_1,
      S_age_2,
      I_age_2,
      S_age_3,
      I_age_3};

/* Offsets in data to handle the infectious pressure */
enum {INFECTIOUS_PRESSURE,
      RESPONSE_age_1,
      RESPONSE_age_2,
      RESPONSE_age_3,
      RECOVER_age_1,
      RECOVER_age_2,
      RECOVER_age_3,
      ALPHA,
      BETA_SEASON_Q1,
      BETA_SEASON_Q2,
      BETA_SEASON_Q3,
      BETA_SEASON_Q4,
      EPSILON};

/**
 * age_1; susceptible to infected: S -> I
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_age_1_to_I_age_1(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return data[RESPONSE_age_1] * data[INFECTIOUS_PRESSURE] * x[S_age_1];
}

/**
 * age_2; susceptible to infected: S -> I
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_age_2_to_I_age_2(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return data[RESPONSE_age_2] * data[INFECTIOUS_PRESSURE] * x[S_age_2];
}

/**
 *  age_3; susceptible to infected: S -> I
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_S_age_3_to_I_age_3(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return data[RESPONSE_age_3] * data[INFECTIOUS_PRESSURE] * x[S_age_3];
}

/**
 *  age_1; infected to susceptible: I -> S
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_I_age_1_to_S_age_1(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return x[I_age_1] / data[RECOVER_age_1];
}

/**
 * age_2; infected to susceptible: I -> S
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe3_I_age_2_to_S_age_2(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return x[I_age_2] / data[RECOVER_age_2];
}

/**
 * age_3; infected to susceptible: I -> S
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity
 */
double SISe3_I_age_3_to_S_age_3(
    const int *x,
    double t,
    const double *data,
    int sd)
{
    return x[I_age_3] / data[RECOVER_age_3];
}

/**
 * Update infectious pressure
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return 1 if needs update, else 0.
 */
int SISe3_update_infectious_pressure(
    const int *x,
    int src,
    double t,
    double *data)
{
    const int days_in_year = 365;
    const int days_in_quarter = 91;

    double S_n, I_n;
    double tmp = data[INFECTIOUS_PRESSURE];

    S_n = x[S_age_1] + x[S_age_2] + x[S_age_3];
    I_n = x[I_age_1] + x[I_age_2] + x[I_age_3];

    /* Time dependent beta */
    data[INFECTIOUS_PRESSURE] *= (1.0 - data[BETA_SEASON_Q1 + ((int)t % days_in_year) / days_in_quarter]);

    data[INFECTIOUS_PRESSURE] += data[ALPHA] * I_n / (I_n + S_n) + data[EPSILON];

    /* 1 if needs update */
    return tmp != data[INFECTIOUS_PRESSURE];
}
