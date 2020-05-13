/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 -- 2020 Stefan Widgren
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

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Visibility.h>

/**
 * Systematic resampling of particles.
 *
 * @param w a numeric vector with weights for the particles.
 * @return an integer vector with indices.
 */
SEXP attribute_hidden SimInf_systematic_resampling(SEXP w)
{
    double cumsum_w, *ptr_w = REAL(w);
    int i, j, n = Rf_length(w);
    double du, u;
    SEXP idx;
    int *ptr_idx;

    for (i = 0, cumsum_w = 0; i < n; i++) {
        if (!R_FINITE(ptr_w[i]) || ptr_w[i] < 0.0)
            Rf_error("Invalid weight detected (non-finite or < 0.0).");
        cumsum_w += ptr_w[i];
    }

    if (cumsum_w <= 0.0)
        Rf_error("Non-positive sum of weights detected.");

    du = cumsum_w / n;
    GetRNGstate();
    u = du * unif_rand();
    PutRNGstate();

    PROTECT(idx = Rf_allocVector(INTSXP, n));
    ptr_idx = INTEGER(idx);
    for (i = 0, j = 0, cumsum_w = ptr_w[0]; i < n; i++) {
        while (u > cumsum_w)
            cumsum_w += ptr_w[++j];
        ptr_idx[i] = j < n ? j + 1 : n;
        u += du;
    }

    UNPROTECT(1);

    return idx;
}
