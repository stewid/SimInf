/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015 - 2017  Stefan Engblom
 *  Copyright (C) 2015 - 2017  Stefan Widgren
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

/**
 * Decay of environmental infectious pressure with a forward Euler step.
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
    double beta_t1, double beta_t2, double beta_t3, double beta_t4)
{
    if (day < end_t2) {
        if (day < end_t1) {
            if (end_t1 < end_t4)
                return phi * (1.0 - beta_t1);
            if (day < end_t4) {
                if (end_t4 < end_t3)
                    return phi * (1.0 - beta_t4);
                if (day < end_t3)
                    return phi * (1.0 - beta_t3);
                return phi * (1.0 - beta_t4);
            }
            return phi * (1.0 - beta_t1);
        }

        return phi * (1.0 - beta_t2);
    }

    if (end_t3 < end_t1 || day < end_t3)
        return phi * (1.0 - beta_t3);

    if (end_t4 < end_t1 || day < end_t4)
        return phi * (1.0 - beta_t4);

    return phi * (1.0 - beta_t1);
}
