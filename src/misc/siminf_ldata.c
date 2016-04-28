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

#include "siminf_arg.h"
#include "siminf_ldata.h"

/**
 * Combine local model parameters and spatial coupling to neighbors in
 * 'ldata'
 *
 * @param data Matrix with local model parameters for each node.
 * @param distance Sparse matrix with distances between nodes.
 * @return matrix
 *
 * Format for ldata:
 * First n indicies (0, 1, ..., n-1): local model parameters followed by
 * the neighbor data, pairs of (index, value) and then a stop pair (-1, 0)
 * where 'index' is the index to the neighbor and value = 1 / distance^2
 *
 * Example:
 * end_t1, end_t2, end_t3, end_t4, index, value, index, value, -1, 0.0
 */
SEXP siminf_ldata_sp(SEXP data, SEXP distance)
{
    SEXP result;
    double *val, *ldata, *ld;
    int i, *ir, *jc, Nld, Nn, node, n_data;

    /* Check arguments */
    if (siminf_arg_check_matrix(data))
        Rf_error("Invalid 'data' argument");
    if (siminf_arg_check_dgCMatrix(distance))
        Rf_error("Invalid 'distance' argument");

    /* Extract data from 'data' */
    Nn = INTEGER(GET_SLOT(data, R_DimSymbol))[1];
    ld = REAL(data);

    /* Extract data from the distance matrix */
    ir = INTEGER(GET_SLOT(distance, Rf_install("i")));
    jc = INTEGER(GET_SLOT(distance, Rf_install("p")));
    val = REAL(GET_SLOT(distance, Rf_install("x")));

    /* Check that the number of nodes are equal in data and
     * distance */
    if (Nn != (LENGTH(GET_SLOT(distance, Rf_install("p"))) - 1))
        Rf_error("The number of nodes in 'data' and 'distance' are not equal");

    /* Calculate length of 'Nld' in 'ldata' for each node in the
     * following three steps: 1), 2), and 3).
     */

    /* 1) Determine the maximum number of neighbors in the 'distance'
     * matrix */
    Nld = 0;
    for (i = 0; i < Nn; i++) {
        const int k = jc[i + 1] - jc[i];
        if (k > Nld)
            Nld = k;
    }

    /* 2) Create one pair for each neighbor (index, value) and add one
     * extra pair for the stop (-1, 0) */
    Nld = (Nld + 1) * 2;

    /*  3) Add space for local model parameters in 'data' */
    n_data = INTEGER(GET_SLOT(data, R_DimSymbol))[0];
    Nld += n_data;

    /* Allocate and initialize memory for ldata */
    PROTECT(result = allocMatrix(REALSXP, Nld, Nn));
    memset(REAL(result), 0, Nn * Nld * sizeof(double));
    ldata = REAL(result);

    for (node = 0; node < Nn; node++) {
        int k = 0;

        /* Copy local model parameters */
        for (i = 0; i < n_data; i++, k++)
            ldata[node * Nld + k] = ld[node * n_data + k];

        /* Copy distance data and recalculate to 1.0 / val^2 */
        for (i = jc[node]; i < jc[node + 1]; i++) {
            ldata[node * Nld + k++] = ir[i];
            ldata[node * Nld + k++] = 1.0 / (val[i] * val[i]);
        }

        /* Add stop */
        ldata[node * Nld + k] = -1.0;
    }

    UNPROTECT(1);

    return result;
}
