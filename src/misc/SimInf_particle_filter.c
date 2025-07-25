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

#include <R.h>
#include <R_ext/Visibility.h>
#include <Rinternals.h>
#include <Rmath.h>

/**
 * Systematic resampling of particles.
 *
 * @param w a numeric vector with weights for the particles.
 * @return an integer vector with indices.
 */
attribute_hidden SEXP SimInf_systematic_resampling(SEXP w)
{
    const double *ptr_w = REAL(w);
    const R_xlen_t n = XLENGTH(w);

    double cumsum_w = 0;
    for (ptrdiff_t i = 0; i < n; i++) {
        if (!R_FINITE(ptr_w[i]) || ptr_w[i] < 0.0)
            Rf_error("Invalid weight detected (non-finite or < 0.0).");
        cumsum_w += ptr_w[i];
    }

    if (cumsum_w <= 0.0)
        Rf_error("Non-positive sum of weights detected.");

    const double du = cumsum_w / (double) n;
    GetRNGstate();
    double u = du * unif_rand();
    PutRNGstate();

    SEXP idx = PROTECT(Rf_allocVector(INTSXP, n));
    int *ptr_idx = INTEGER(idx);
    cumsum_w = ptr_w[0];
    for (ptrdiff_t i = 0, j = 0; i < n; i++) {
        while (u > cumsum_w)
            cumsum_w += ptr_w[++j];
        ptr_idx[i] = (int) (j < n ? j + 1 : n);
        u += du;
    }

    UNPROTECT(1);

    return idx;
}

/**
 * Split scheduled events into the intervals in tspan.
 *
 * @param t an integer vector with the time-points for the events.
 * @param t_end an integer vector with the time end-point in each
 *        interval.
 * @return a len(t_end) x 2 integer matrix where the first column
 *         contains the one-based index to the first event in the
 *         interval, and the second column the number of events in the
 *         interval.
 */
attribute_hidden SEXP SimInf_split_events(SEXP t, SEXP t_end)
{
    const int *ptr_t = INTEGER(t);
    const R_xlen_t t_len = XLENGTH(t);
    if (t_len < 1)
        Rf_error("'t' must be an integer vector with length >= 1.");

    const int *ptr_t_end = INTEGER(t_end);
    const R_xlen_t t_end_len = XLENGTH(t_end);
    if (t_end_len < 1)
        Rf_error("'t_end' must be an integer vector with length >= 1.");

    /* Create a matrix to hold the result. */
    SEXP m = PROTECT(Rf_allocMatrix(INTSXP, t_end_len, 2));
    int *ptr_m = INTEGER(m);
    memset(ptr_m, 0, t_end_len * 2 * sizeof(int));

    /* Interate over the event time-points and place each event in the
     * corresponding interval. */
    ptrdiff_t t_i = 0;
    ptrdiff_t t_end_i = 0;
    while (t_i < t_len && t_end_i < t_end_len) {
        if (ptr_t[t_i] <= ptr_t_end[t_end_i]) {
            if (!ptr_m[t_end_i])
                ptr_m[t_end_i] = (int) (t_i + 1);
            ptr_m[t_end_len + t_end_i]++;
            t_i++;
        } else {
            t_end_i++;
        }
    }

    UNPROTECT(1);

    return m;
}
