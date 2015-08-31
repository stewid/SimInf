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

#include "siminf_events.h"

/**
 * Create a directed random graph G_{n,p} without loops
 *
 * Reference: Batagelj, V. & Brandes, U. Efficient generation of large
 * random networks Physical Review E, APS, 2005, 71, 036113
 * @param n Number of nodes
 * @param p Probability of edge
 * @param from The vector to add the 'from' nodes to
 * @param to The vector to add the 'to' nodes to
 * @param rng The random number generator
 * @return 0 if Ok, else error code.
*/
static int siminf_gnp(
    int n,
    double p,
    struct siminf_vec *from,
    struct siminf_vec *to,
    gsl_rng *rng)
{
    int w = -1;
    int v = 0;
    double lp;

    if (p <= 0.0 || p >= 1.0)
        return SIMINF_INVALID_EDGE_PROBABILITY;

    lp = log(1.0 - p);

    while (v < n) {
        double lr = log(1.0 - gsl_rng_uniform(rng));
        w = w + 1 + lr / lp;

        /* Protect against loops */
        if (v == w)
            w += 1;

        while  (w >= n && v < n) {
            w = w - v;
            v = v + 1;

            /* Protect against loops */
            if (v == w)
                w += 1;
        }

        if (v < n) {
            int err;

            err = siminf_vec_push_back(from, v);
            if (err)
                return err;

            err = siminf_vec_push_back(to, w);
            if (err)
                return err;
        }
    }

    return 0;
}
