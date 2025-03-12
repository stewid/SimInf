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

#include <Rinternals.h>
#include <R_ext/Visibility.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_lambert.h>

/**
 * Utility function to compute the principal branch of the Lambert W
 * function, W_0(x).
 *
 * @param x a numeric vector.
 * @return a numeric vector with the same length as the input vector x.
 */
attribute_hidden
SEXP
SimInf_lambertW0(
    SEXP x)
{
    SEXP W0;
    R_xlen_t len;

    if (!Rf_isReal(x))
        Rf_error("'x' must be a numeric vector.");

    len = XLENGTH(x);
    PROTECT(W0 = Rf_allocVector(REALSXP, len));

    for (R_xlen_t i = 0; i < len; i++) {
        double xx= REAL(x)[i];
        double val = R_NaN;
        gsl_sf_result result;

        if (ISNA(xx))
            val = NA_REAL;
        else if (xx == R_PosInf)
            val = R_PosInf;
        else if (R_FINITE(xx) && gsl_sf_lambert_W0_e(xx, &result) == GSL_SUCCESS)
            val = result.val;
        REAL(W0)[i] = val;
    }

    UNPROTECT(1);

    return W0;
}
