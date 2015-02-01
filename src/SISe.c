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
enum {S, I};

/* Offsets in data to handle the infectious pressure */
enum {INFECTIOUS_PRESSURE,
      RESPONSE,
      RECOVER,
      ALPHA,
      BETA_SEASON_Q1,
      BETA_SEASON_Q2,
      BETA_SEASON_Q3,
      BETA_SEASON_Q4,
      EPSILON};

/**
 * susceptible to infected: S -> I
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe_S_to_I(const int *x, double t, const double *data, int sd)
{
    return data[RESPONSE] * data[INFECTIOUS_PRESSURE] * x[S];
}

/**
 *  infected to susceptible: I -> S
 *
 * @param x The state vector in node.
 * @param t Current time.
 * @param data The data vector for node.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe_I_to_S(const int *x, double t, const double *data, int sd)
{
    return x[I] / data[RECOVER];
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
int SISe_update_infectious_pressure(
    const int *x,
    int src,
    double t,
    double *data)
{
    const int days_in_year = 365;
    const int days_in_quarter = 91;

    double S_n, I_n;
    double tmp = data[INFECTIOUS_PRESSURE];

    S_n = x[S];
    I_n = x[I];

    /* Time dependent beta */
    data[INFECTIOUS_PRESSURE] *= (1.0 - data[BETA_SEASON_Q1 + ((int)t % days_in_year) / days_in_quarter]);

    data[INFECTIOUS_PRESSURE] += data[ALPHA] * I_n / (I_n + S_n) + data[EPSILON];

    /* 1 if needs update */
    return tmp != data[INFECTIOUS_PRESSURE];
}
