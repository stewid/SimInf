/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
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

#include "SimInf.h"
#include <R_ext/Visibility.h>
#include <Rinternals.h>

static ptrdiff_t
SimInf_Euclidean_distance(const double *x,
                          const double *y,
                          const double cutoff,
                          const double min_dist,
                          const int na_fail,
                          const ptrdiff_t len,
                          double *distance,
                          int *row_indices, int *col_indices)
{
    if (col_indices)
        col_indices[0] = 0;

    ptrdiff_t n = 0;
    for (ptrdiff_t i = 0; i < len; i++) {
        for (ptrdiff_t j = 0; j < len; j++) {
            if (i != j) {
                /* Calculate the Euclidean distance. */
                double d = hypot(x[i] - x[j], y[i] - y[j]);

                if (!R_FINITE(d)) {
                    if ((R_IsNA(x[i]) ||
                         R_IsNA(x[j]) ||
                         R_IsNA(y[i]) || R_IsNA(y[j])) && !na_fail) {
                        continue;
                    }

                    Rf_error("Invalid distance for i=%"
                             R_PRIdXLEN_T
                             " and j=%" R_PRIdXLEN_T ".", i, j);
                }

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
                        row_indices[n] = (int) j;

                    n++;
                }
            }
        }

        if (col_indices)
            col_indices[i + 1] = (int) n;
    }

    return n;
}

attribute_hidden
    SEXP
SimInf_distance_matrix(SEXP x_,
                       SEXP y_,
                       SEXP cutoff_, SEXP min_dist_, SEXP na_fail_)
{
    /* Check that the input vectors have an identical length > 0. */
    const R_xlen_t len = XLENGTH(x_);
    if (len < 1)
        Rf_error("'x' must be a numeric vector with length > 0.");
    if (XLENGTH(y_) != len) {
        Rf_error("'y' must be a numeric vector with length %" R_PRIdXLEN_T
                 ".", len);
    }

    /* Check for valid cutoff. */
    const double cutoff = Rf_asReal(cutoff_);
    if (!R_FINITE(cutoff) || cutoff < 0)
        Rf_error("'cutoff' must be > 0.");

    /* Check for a valid na_fail. */
    if (!Rf_isLogical(na_fail_) ||
        Rf_length(na_fail_) != 1 || LOGICAL(na_fail_)[0] == NA_LOGICAL) {
        Rf_error("'na_fail' must be TRUE or FALSE.");
    }
    const int na_fail = LOGICAL(na_fail_)[0];

    /* First, iterate over all the elements to determine the required
     * length for the result vector. */
    const double *x = REAL(x_);
    const double *y = REAL(y_);
    const double min_dist = Rf_asReal(min_dist_);
    const ptrdiff_t n = SimInf_Euclidean_distance(x,
                                                  y,
                                                  cutoff,
                                                  min_dist,
                                                  na_fail,
                                                  len,
                                                  NULL,
                                                  NULL,
                                                  NULL);

    /* Allocate vectors for the sparse matrix. */
    SEXP distance = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP row_indices = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP col_indices = PROTECT(Rf_allocVector(INTSXP, len + 1));

    /* Now, iterate over all the elements again and save the result in
     * the allocated result vectors. */
    SimInf_Euclidean_distance(x,
                              y,
                              cutoff,
                              min_dist,
                              na_fail,
                              len,
                              REAL(distance),
                              INTEGER(row_indices), INTEGER(col_indices));

    /* Create the sparse matrix. */
    SEXP class = PROTECT(R_do_MAKE_CLASS("dgCMatrix"));
    SEXP result = PROTECT(R_do_new_object(class));
    R_do_slot_assign(result, Rf_install("x"), distance);
    R_do_slot_assign(result, Rf_install("i"), row_indices);
    R_do_slot_assign(result, Rf_install("p"), col_indices);
    INTEGER(R_do_slot(result, Rf_install("Dim")))[0] = (int) len;
    INTEGER(R_do_slot(result, Rf_install("Dim")))[1] = (int) len;

    UNPROTECT(5);

    return result;
}
