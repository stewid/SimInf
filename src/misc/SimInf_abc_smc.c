/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 -- 2019 Stefan Widgren
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "SimInf_arg.h"

/**
 * Utility function for implementing the Approximate Bayesian
 * Computation Sequential Monte Carlo (ABC-SMC) algorithm of Toni et
 * al. (2009). Samples particles from: first the prior distribution
 * for the first generation, and then from a previous generation of
 * particles.
 *
 * @param parameter character vector with the name of the parameter
 *        for each prior.
 * @param distribution character vector with the name of the
 *        distribution for each prior. Each entry must contain one of
 *        'G' (gamma), 'N' (normal) or 'U' (uniform).
 * @param p1 numeric vector with the first hyperparameter for each
 *        prior: G) shape, N) mean, and U) lower bound.
 * @param p2 numeric vector with the second hyperparameter for each
 *        prior: G) rate, N) standard deviation, and U) upper bound.
 * @param n number of proposals to generate.
 * @param x a numeric matrix (parameters x particles) with a previous
 *        generation of particles or NULL.
 * @param w a numeric vector with weigths for the previous generation
 *        of particles or NULL.
 * @param sigma variance-covariance matrix.
 * @return a numeric matrix (parameters x particles) with
 *         proposals. The matrix also has an attribute 'ancestor' with
 *         an index that indicates which particle it was sampled from.
 */
SEXP SimInf_abc_smc_proposals(
    SEXP parameter,
    SEXP distribution,
    SEXP p1,
    SEXP p2,
    SEXP n,
    SEXP x,
    SEXP w,
    SEXP sigma)
{
    int error = 0, k, len, N;
    gsl_rng *rng = NULL;
    gsl_matrix_view v_sigma;
    gsl_matrix *SIGMA = NULL;
    double *ptr_x = NULL, *ptr_w = NULL, *cdf = NULL;
    double *ptr_p1 = REAL(p1), *ptr_p2 = REAL(p2);
    SEXP xx, ancestor, dimnames;
    double *ptr_xx;
    int *ptr_ancestor;

    /* Check input arguments. */
    if (SimInf_arg_check_integer_gt_zero(n))
        Rf_error("'n' must be an integer > 0.");
    N = INTEGER(n)[0];
    if (!Rf_isString(parameter))
        Rf_error("'parameter' must be a character vector.");
    k = Rf_length(parameter);

    /* Setup result matrix. */
    PROTECT(xx = Rf_allocMatrix(REALSXP, k, N));
    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    Rf_setAttrib(xx, R_DimNamesSymbol, dimnames);
    SET_VECTOR_ELT(dimnames, 0, parameter);
    ptr_xx = REAL(xx);

    /* Setup vector to record 'ancestor' i.e. which particle the
     * proposal was sampled from. */
    PROTECT(ancestor = Rf_allocVector(INTSXP, N));
    Rf_setAttrib(xx, Rf_install("ancestor"), ancestor);
    ptr_ancestor = INTEGER(ancestor);

    /* Setup random number generator. */
    GetRNGstate();
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!rng) {
        error = 1;
        goto cleanup;
    }
    gsl_rng_set(rng, runif(0.0, 1.0) * UINT_MAX);

    if (Rf_isNull(x)) {
        /* First generation: sample from priors. */
        for (int i = 0; i < N; i++) {
            ptr_ancestor[i] = NA_INTEGER;
            for (int d = 0; d < k; d++) {
                switch(R_CHAR(STRING_ELT(distribution, d))[0]) {
                case 'G':
                    ptr_xx[i * k + d] = rgamma(ptr_p1[d], 1.0 / ptr_p2[d]);
                    break;
                case 'N':
                    ptr_xx[i * k + d] = rnorm(ptr_p1[d], ptr_p2[d]);
                    break;
                case 'U':
                    ptr_xx[i * k + d] = runif(ptr_p1[d], ptr_p2[d]);
                    break;
                default:
                    error = 2;
                    break;
                }
            }
        }

        goto cleanup;
    }

    /* Setup variance-covariance matrix. */
    v_sigma = gsl_matrix_view_array(REAL(sigma), k, k);
    SIGMA = gsl_matrix_alloc(k, k);
    if (!SIGMA) {
        error = 1;
        goto cleanup;
    }
    gsl_matrix_memcpy(SIGMA, &v_sigma.matrix);
    gsl_linalg_cholesky_decomp1(SIGMA);

    /* Setup weights */
    ptr_x = REAL(x);
    ptr_w = REAL(w);
    len = Rf_length(w);
    cdf = malloc(len * sizeof(double));
    for (int i = 0; i < len; i++) {
        if (!R_FINITE(ptr_w[i]) || ptr_w[i] < 0.0) {
            error = 3;
            goto cleanup;
        }
        cdf[i] = ptr_w[i];
        if (i > 0)
            cdf[i] += cdf[i-1];
    }

    for (int i = 0; i < N; i++) {
        int accept;
        gsl_vector_view X;
        gsl_vector_view proposal = gsl_vector_view_array(&ptr_xx[i * k], k);

        do {
            /* Sample a particle from previous generation. Use a
             * binary search to determine the sampled particle based
             * on its weight: [0, cdf_0), [cdf_0, cdf_1), ... */
            int j = 0, j_low = 0, j_high = len - 1;
            double r = unif_rand() * cdf[j_high]; /* r ~ U[0, cdf[j_high]) */
            while (j_low < j_high) {
                j = (j_low + j_high) / 2;
                if (cdf[j] <= r)
                    j_low = j + 1;
                else
                    j_high = j - 1;
            }
            ptr_ancestor[i] = j + 1; /* R is one-based. */

            /* Perturbate the particle. */
            X = gsl_vector_view_array(&ptr_x[j * k], k);
            gsl_ran_multivariate_gaussian(rng, &X.vector, SIGMA,
                                          &proposal.vector);

            /* Check that the proposal is valid. */
            accept = 1;
            for (int d = 0; d < k; d++) {
                switch(R_CHAR(STRING_ELT(distribution, d))[0]) {
                case 'G':
                    if (dgamma(ptr_xx[i * k + d], ptr_p1[d],
                               1.0 / ptr_p2[d], 0) == 0) {
                        accept = 0;
                    }
                    break;
                case 'N':
                    if (dnorm(ptr_xx[i * k + d], ptr_x[j * k + d],
                              ptr_p2[d], 0) == 0) {
                        accept = 0;
                    }
                    break;
                case 'U':
                    if (dunif(ptr_xx[i * k + d], ptr_p1[d],
                              ptr_p2[d], 0) == 0) {
                        accept = 0;
                    }
                    break;
                default:
                    error = 2;
                    break;
                }
            }
        } while(!accept);
    }

cleanup:
    free(cdf);
    gsl_matrix_free(SIGMA);
    gsl_rng_free(rng);
    PutRNGstate();

    if (error) {
        switch (error) {
        case 1:
            Rf_error("Unable to allocate memory buffer.");
            break;
        case 2:
            Rf_error("Unknown distribution.");
            break;
        case 3:
            Rf_error("Invalid weight detected (non-finite or < 0.0).");
            break;
        default:
            Rf_error("Unknown error code: %i.", error);
            break;
        }
    }

    UNPROTECT(3);

    return xx;
}

/**
 * Utility function for implementing the Approximate Bayesian
 * Computation Sequential Monte Carlo (ABC-SMC) algorithm of Toni et
 * al. (2009). Calculate weights for current generation of particles.
 *
 * @param distribution character vector with the name of the
 *        distribution for each prior. Each entry must contain one of
 *        'G' (gamma), 'N' (normal) or 'U' (uniform).
 * @param p1 numeric vector with the first hyperparameter for each
 *        prior: G) shape, N) mean, and U) lower bound.
 * @param p2 numeric vector with the second hyperparameter for each
 *        prior: G) rate, N) standard deviation, and U) upper bound.
 * @param x a numeric matrix (parameters x particles) with the
 *        previous generation of particles or NULL.
 * @param xx a numeric matrix (parameters x particles) with the
 *        current generation of particles or NULL.
 * @param w a numeric vector with weights for the previous generation
 *        of particles.
 * @param sigma variance-covariance matrix.
 * @return a numeric vector with weights for the current generation of
 *         particles.
 */
SEXP SimInf_abc_smc_weights(
    SEXP distribution,
    SEXP p1,
    SEXP p2,
    SEXP x,
    SEXP xx,
    SEXP w,
    SEXP sigma)
{
    int error = 0;
    int k, n = Rf_ncols(xx);
    gsl_matrix_view v_sigma;
    gsl_matrix *SIGMA = NULL;
    gsl_vector *work = NULL;
    SEXP ww;
    double *ptr_p1, *ptr_p2, *ptr_x, *ptr_xx, *ptr_w, *ptr_ww;
    double sum = 0.0, max_ww = 0.0;

    PROTECT(ww = Rf_allocVector(REALSXP, n));
    ptr_ww = REAL(ww);
    if (Rf_isNull(w)) {
        for (int i = 0; i < n; ++i)
            ptr_ww[i] = 1.0 / (double)n;
        goto cleanup;
    }

    k = INTEGER(GET_SLOT(sigma, R_DimSymbol))[0];
    ptr_p1 = REAL(p1);
    ptr_p2 = REAL(p2);
    ptr_x = REAL(x);
    ptr_xx = REAL(xx);
    ptr_w = REAL(w);
    work = gsl_vector_alloc(k);

    /* Setup variance-covariance matrix. */
    v_sigma = gsl_matrix_view_array(REAL(sigma), k, k);
    SIGMA = gsl_matrix_alloc(k, k);
    if (!SIGMA) {
        error = 1;
        goto cleanup;
    }
    gsl_matrix_memcpy(SIGMA, &v_sigma.matrix);
    gsl_linalg_cholesky_decomp1(SIGMA);

    for (int i = 0; i < n; i++) {
        gsl_vector_view v_x = gsl_vector_view_array(&ptr_x[i * k], k);
        gsl_vector_view v_xx = gsl_vector_view_array(&ptr_xx[i * k], k);
        double pdf;

        ptr_ww[i] = 0.0;
        for (int d = 0; d < k; d++) {
            switch(R_CHAR(STRING_ELT(distribution, d))[0]) {
            case 'G':
                ptr_ww[i] +=
                    dgamma(ptr_xx[i * k + d], ptr_p1[d], 1.0 / ptr_p2[d], 1);
                break;
            case 'N':
                ptr_ww[i] +=
                    dnorm(ptr_xx[i * k + d], ptr_x[i * k + d], ptr_p2[d], 1);
                break;
            case 'U':
                ptr_ww[i] +=
                    dunif(ptr_xx[i * k + d], ptr_p1[d], ptr_p2[d], 1);
                break;
            default:
                error = 2;
                break;
            }
        }

        gsl_ran_multivariate_gaussian_pdf(&v_xx.vector, &v_x.vector,
                                          SIGMA, &pdf, work);
        if (!R_FINITE(ptr_w[i]) || ptr_w[i] < 0.0) {
            error = 3;
            goto cleanup;
        }
        sum += ptr_w[i] * pdf;
    }

    sum = log(sum);
    for (int i = 0; i < n; i++) {
        ptr_ww[i] -= sum;
        if (ptr_ww[i] > max_ww)
            max_ww = ptr_ww[i];
    }

    sum = 0.0;
    for (int i = 0; i < n; i++) {
        ptr_ww[i] -= max_ww;
        ptr_ww[i] = exp(ptr_ww[i]);
        sum += ptr_ww[i];
    }

    /* Normalize. */
    for (int i = 0; i < n; i++)
        ptr_ww[i] /= sum;

cleanup:
    gsl_matrix_free(SIGMA);
    gsl_vector_free(work);

    if (error) {
        switch (error) {
        case 1:
            Rf_error("Unable to allocate memory buffer.");
            break;
        case 2:
            Rf_error("Unknown distribution.");
            break;
        case 3:
            Rf_error("Invalid weight detected (non-finite or < 0.0).");
            break;
        default:
            Rf_error("Unknown error code: %i.", error);
            break;
        }
    }

    UNPROTECT(1);

    return ww;
}
