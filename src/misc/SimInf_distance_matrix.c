/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 -- 2022 Stefan Widgren
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

#include <Rdefines.h>
#include <R_ext/Visibility.h>

static R_xlen_t
SimInf_Euclidean_distance(
    const double* x,
    const double* y,
    double cutoff,
    double min_dist,
    R_xlen_t len,
    double *distance,
    int *row_indices,
    int *col_indices)
{
    R_xlen_t n = 0;

    if (col_indices)
        col_indices[0] = 0;

    for (R_xlen_t i = 0; i < len; i++) {
        for (R_xlen_t j = 0; j < len; j++) {
            if (i != j) {
                /* Calculate the Euclidean distance. */
                double d = hypot(x[i] - x[j], y[i] - y[j]);

                if (!R_FINITE(d))
                    Rf_error("Invalid distance for i=%i and j=%i.", i, j);

                if (d <= cutoff) {
                    if (d <= 0) {
                        if (!R_FINITE(min_dist) || min_dist < 0) {
                            Rf_error("Invalid 'min_dist' argument. "
                                     "Please provide 'min_dist' > 0.");
                        }

                        d = min_dist;
                    }

                    if (distance)
                        distance[n] = d;

                    if (row_indices)
                        row_indices[n] = j;

                    n++;
                }
            }
        }

        if (col_indices)
            col_indices[i + 1] = n;
    }

    return n;
}

SEXP attribute_hidden
SimInf_distance_matrix(
    SEXP x_,
    SEXP y_,
    SEXP cutoff_,
    SEXP min_dist_)
{
    double *x = REAL(x_);
    double *y = REAL(y_);
    double cutoff = Rf_asReal(cutoff_);
    double min_dist = Rf_asReal(min_dist_);
    R_xlen_t len = XLENGTH(x_);
    R_xlen_t n;
    SEXP distance;
    SEXP row_indices;
    SEXP col_indices;
    SEXP result;

    /* Check that the input vectors have an identical length > 0. */
    if (len < 1)
        Rf_error("'x' must be a numeric vector with length > 0.");
    if (XLENGTH(y_) != len)
        Rf_error("'y' must be a numeric vector with length %i.", len);

    /* Check for valid cutoff. */
    if (!R_FINITE(cutoff) || cutoff < 0)
        Rf_error("'cutoff' must be > 0.");

    /* First, iterate over all the elements to determine the required
     * length for the result vector. */
    n = SimInf_Euclidean_distance(
        x,
        y,
        cutoff,
        min_dist,
        len,
        NULL,
        NULL,
        NULL);

    /* Allocate vectors for the sparse matrix. */
    PROTECT(distance = Rf_allocVector(REALSXP, n));
    PROTECT(row_indices = Rf_allocVector(INTSXP, n));
    PROTECT(col_indices = Rf_allocVector(INTSXP, len + 1));

    /* Now, iterate over all the elements again and save the result in
     * the allocated result vectors. */
    SimInf_Euclidean_distance(
        x,
        y,
        cutoff,
        min_dist,
        len,
        REAL(distance),
        INTEGER(row_indices),
        INTEGER(col_indices));

    /* Create the sparse matrix. */
    PROTECT(result = NEW_OBJECT(MAKE_CLASS("dgCMatrix")));
    SET_SLOT(result, Rf_install("x"), distance);
    SET_SLOT(result, Rf_install("i"), row_indices);
    SET_SLOT(result, Rf_install("p"), col_indices);
    INTEGER(GET_SLOT(result, Rf_install("Dim")))[0] = len;
    INTEGER(GET_SLOT(result, Rf_install("Dim")))[1] = len;

    UNPROTECT(4);

    return result;
}
