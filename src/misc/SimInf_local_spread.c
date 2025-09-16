/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
 * Copyright (C) 2015 -- 2025 Stefan Widgren
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

#include <R_ext/Visibility.h>

/**
 * Local spread of the environmental infectious pressure phi among
 * proximal nodes.
 *
 * @param neighbors Spatial coupling between nodes where 'neighbors'
 * is a vector of pairs (index, distance) to neighbor nodes. The pair
 * vector is terminated with an index equal to -1.
 * @param phi Vector with phi in each node
 * @param u The compartment state vector in each node.
 * @param N_i The number of individuals in node i.
 * @param phi_i The environmental infectious pressure phi in node i.
 * @param Nc The number of compartments in each node.
 * @param D The spatial coupling of the environmental infectious
 * pressure phi among proximal nodes.
 * @return The contribution from neighbors to phi in node i
 */
attribute_hidden double
SimInf_local_spread(
    const double *neighbors,
    const double *phi,
    const int *u,
    const double N_i,
    const double phi_i,
    const int Nc,
    const double D)
{
    double ls = 0.0;
    const double phi_i_N_i = phi_i * N_i;

    int j = (int) *neighbors++;
    while (j >= 0) {
        /* Count number of individuals in node j */
        double N_j = 0;
        for (int k = j * Nc; k < (j + 1) * Nc; k++)
            N_j += u[k];

        if (N_j > 0.0)
            ls += ((phi[j] * N_j - phi_i_N_i) * D) / (N_i * (*neighbors));

        /* Move to next neighbor pair (index, distance) */
        neighbors++;
        j = (int) *neighbors++;
    }

    return ls;
}
