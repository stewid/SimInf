/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 -- 2021 Stefan Widgren
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

#if defined SIMINF_NO_ABC

# include <Rdefines.h>
# include <R_ext/Visibility.h>
# include "SimInf.h"

SEXP attribute_hidden SimInf_abc_proposals(
    SEXP parameter,
    SEXP distribution,
    SEXP p1,
    SEXP p2,
    SEXP n,
    SEXP x,
    SEXP w,
    SEXP sigma)
{
    SIMINF_UNUSED(parameter);
    SIMINF_UNUSED(distribution);
    SIMINF_UNUSED(p1);
    SIMINF_UNUSED(p2);
    SIMINF_UNUSED(n);
    SIMINF_UNUSED(x);
    SIMINF_UNUSED(w);
    SIMINF_UNUSED(sigma);

    Rf_error("The installed version of the GNU Scientific Library (GSL) that "
             "is required to build SimInf with support for ABC is to old. "
             "Please install GSL version >= 2.2 and reinstall SimInf if you "
             "need that functionality.");

    return R_NilValue;
}

SEXP attribute_hidden SimInf_abc_weights(
    SEXP distribution,
    SEXP p1,
    SEXP p2,
    SEXP x,
    SEXP xx,
    SEXP w,
    SEXP sigma)
{
    SIMINF_UNUSED(distribution);
    SIMINF_UNUSED(p1);
    SIMINF_UNUSED(p2);
    SIMINF_UNUSED(x);
    SIMINF_UNUSED(xx);
    SIMINF_UNUSED(w);
    SIMINF_UNUSED(sigma);

    Rf_error("The installed version of the GNU Scientific Library (GSL) that "
             "is required to build SimInf with support for ABC is to old. "
             "Please install GSL version >= 2.2 and reinstall SimInf if you "
             "need that functionality.");

    return R_NilValue;
}

#else

# include <R.h>
# include <Rdefines.h>
# include <Rmath.h>
# include <R_ext/Visibility.h>
# include <gsl/gsl_matrix.h>
# include <gsl/gsl_linalg.h>
# include <gsl/gsl_randist.h>
# include <gsl/gsl_rng.h>
# include "SimInf_arg.h"

static void SimInf_abc_error(int error)
{
    switch (error) {
    case 1:                                            /* #nocov */
        Rf_error("Unable to allocate memory buffer."); /* #nocov */
        break;
    case 2:
        Rf_error("Unknown distribution.");
        break;
    case 3:
        Rf_error("Invalid weight detected (non-finite or < 0.0).");
        break;
    default:                                        /* #nocov */
        Rf_error("Unknown error code: %i.", error); /* #nocov */
        break;
    }
}

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
 *        'g' (gamma), 'n' (normal) or 'u' (uniform).
 * @param p1 numeric vector with the first hyperparameter for each
 *        prior: g) shape, n) mean, and u) lower bound.
 * @param p2 numeric vector with the second hyperparameter for each
 *        prior: g) rate, n) standard deviation, and u) upper bound.
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
SEXP attribute_hidden SimInf_abc_proposals(
    SEXP parameter,
    SEXP distribution,
    SEXP p1,
    SEXP p2,
    SEXP n,
    SEXP x,
    SEXP w,
    SEXP sigma)
{
    int error = 0, n_parameters, len = 0, n_proposals;
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
    n_proposals = INTEGER(n)[0];
    if (!Rf_isString(parameter))
        Rf_error("'parameter' must be a character vector.");
    n_parameters = Rf_length(parameter);
    if (!Rf_isNull(x)) {
        len = Rf_length(w);
        if (len < 1)
            Rf_error("'w' must have length >= 1 when 'x' is non-null.");
    }

    /* Setup result matrix. */
    PROTECT(xx = Rf_allocMatrix(REALSXP, n_parameters, n_proposals));
    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    Rf_setAttrib(xx, R_DimNamesSymbol, dimnames);
    SET_VECTOR_ELT(dimnames, 0, parameter);
    ptr_xx = REAL(xx);

    /* Setup vector to record 'ancestor' i.e. which particle the
     * proposal was sampled from. */
    PROTECT(ancestor = Rf_allocVector(INTSXP, n_proposals));
    Rf_setAttrib(xx, Rf_install("ancestor"), ancestor);
    ptr_ancestor = INTEGER(ancestor);

    /* Setup random number generator. */
    GetRNGstate();
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    if (!rng) {
        error = 1;    /* #nocov */
        goto cleanup; /* #nocov */
    }
    gsl_rng_set(rng, runif(1, UINT_MAX));

    if (Rf_isNull(x)) {
        /* First generation: sample from priors. */
        for (int i = 0; i < n_proposals; i++) {
            ptr_ancestor[i] = NA_INTEGER;
            for (int d = 0; d < n_parameters; d++) {
                switch(R_CHAR(STRING_ELT(distribution, d))[0]) {
                case 'g':
                    ptr_xx[i * n_parameters + d] =
                        rgamma(ptr_p1[d], 1.0 / ptr_p2[d]);
                    break;
                case 'n':
                    ptr_xx[i * n_parameters + d] = rnorm(ptr_p1[d], ptr_p2[d]);
                    break;
                case 'u':
                    ptr_xx[i * n_parameters + d] = runif(ptr_p1[d], ptr_p2[d]);
                    break;
                default:
                    error = 2;
                    goto cleanup;
                }
            }
        }

        goto cleanup;
    }

    /* Setup variance-covariance matrix. */
    v_sigma = gsl_matrix_view_array(REAL(sigma), n_parameters, n_parameters);
    SIGMA = gsl_matrix_alloc(n_parameters, n_parameters);
    if (!SIGMA) {
        error = 1;    /* #nocov */
        goto cleanup; /* #nocov */
    }
    gsl_matrix_memcpy(SIGMA, &v_sigma.matrix);
    gsl_linalg_cholesky_decomp1(SIGMA);

    /* Setup weights */
    ptr_x = REAL(x);
    ptr_w = REAL(w);
    cdf = malloc(len * sizeof(double));
    if (!cdf) {
        error = 1;    /* #nocov */
        goto cleanup; /* #nocov */
    }
    for (int i = 0; i < len; i++) {
        if (!R_FINITE(ptr_w[i]) || ptr_w[i] < 0.0) {
            error = 3;
            goto cleanup;
        }
        cdf[i] = ptr_w[i];
        if (i > 0)
            cdf[i] += cdf[i-1];
    }

    for (int i = 0; i < n_proposals; i++) {
        int accept;
        gsl_vector_view X;
        gsl_vector_view proposal =
            gsl_vector_view_array(&ptr_xx[i * n_parameters], n_parameters);

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
            X = gsl_vector_view_array(&ptr_x[j * n_parameters], n_parameters);
            gsl_ran_multivariate_gaussian(rng, &X.vector, SIGMA,
                                          &proposal.vector);

            /* Check that the proposal is valid. */
            accept = 1;
            for (int d = 0; d < n_parameters; d++) {
                double density;

                switch(R_CHAR(STRING_ELT(distribution, d))[0]) {
                case 'g':
                    density = dgamma(ptr_xx[i * n_parameters + d],
                                     ptr_p1[d],
                                     1.0 / ptr_p2[d],
                                     0);
                    break;
                case 'n':
                    density = dnorm(ptr_xx[i * n_parameters + d],
                                    ptr_x[j * n_parameters + d],
                                    ptr_p2[d],
                                    0);
                    break;
                case 'u':
                    density = dunif(ptr_xx[i * n_parameters + d],
                                    ptr_p1[d],
                                    ptr_p2[d],
                                    0);
                    break;
                default:
                    error = 2;
                    goto cleanup;
                }

                if (!R_FINITE(density) || density <= 0.0)
                    accept = 0;
            }
        } while(!accept);
    }

cleanup:
    free(cdf);
    gsl_matrix_free(SIGMA);
    gsl_rng_free(rng);
    PutRNGstate();

    if (error)
        SimInf_abc_error(error);

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
 *        'g' (gamma), 'n' (normal) or 'u' (uniform).
 * @param p1 numeric vector with the first hyperparameter for each
 *        prior: g) shape, n) mean, and u) lower bound.
 * @param p2 numeric vector with the second hyperparameter for each
 *        prior: g) rate, n) standard deviation, and u) upper bound.
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
SEXP attribute_hidden SimInf_abc_weights(
    SEXP distribution,
    SEXP p1,
    SEXP p2,
    SEXP x,
    SEXP xx,
    SEXP w,
    SEXP sigma)
{
    int error = 0;
    int n_parameters, n_particles = Rf_ncols(xx);
    gsl_matrix_view v_sigma;
    gsl_matrix *SIGMA = NULL;
    gsl_vector *work = NULL;
    SEXP ww;
    double *ptr_p1, *ptr_p2, *ptr_x, *ptr_xx, *ptr_w, *ptr_ww;
    double sum, max_ww = 0.0;

    PROTECT(ww = Rf_allocVector(REALSXP, n_particles));
    ptr_ww = REAL(ww);
    if (Rf_isNull(w)) {
        for (int i = 0; i < n_particles; ++i)
            ptr_ww[i] = 1.0 / (double)n_particles;
        goto cleanup;
    }

    n_parameters = INTEGER(GET_SLOT(sigma, R_DimSymbol))[0];
    ptr_p1 = REAL(p1);
    ptr_p2 = REAL(p2);
    ptr_x = REAL(x);
    ptr_xx = REAL(xx);
    ptr_w = REAL(w);
    work = gsl_vector_alloc(n_parameters);

    /* Setup variance-covariance matrix. */
    v_sigma = gsl_matrix_view_array(REAL(sigma), n_parameters, n_parameters);
    SIGMA = gsl_matrix_alloc(n_parameters, n_parameters);
    if (!SIGMA) {
        error = 1;    /* #nocov */
        goto cleanup; /* #nocov */
    }
    gsl_matrix_memcpy(SIGMA, &v_sigma.matrix);
    gsl_linalg_cholesky_decomp1(SIGMA);

    for (int i = 0; i < n_particles; i++) {
        gsl_vector_view v_xx = gsl_vector_view_array(&ptr_xx[i * n_parameters],
                                                     n_parameters);

        ptr_ww[i] = 0.0;
        for (int d = 0; d < n_parameters; d++) {
            switch(R_CHAR(STRING_ELT(distribution, d))[0]) {
            case 'g':
                ptr_ww[i] +=
                    dgamma(ptr_xx[i * n_parameters + d],
                           ptr_p1[d], 1.0 / ptr_p2[d], 1);
                break;
            case 'n':
                ptr_ww[i] +=
                    dnorm(ptr_xx[i * n_parameters + d],
                          ptr_x[i * n_parameters + d], ptr_p2[d], 1);
                break;
            case 'u':
                ptr_ww[i] +=
                    dunif(ptr_xx[i * n_parameters + d],
                          ptr_p1[d], ptr_p2[d], 1);
                break;
            default:
                error = 2;
                goto cleanup;
            }
        }

        if (!R_FINITE(ptr_w[i]) || ptr_w[i] < 0.0) {
            error = 3;
            goto cleanup;
        }

        sum = 0.0;
        for (int j = 0; j < n_particles; j++) {
            double pdf;
            gsl_vector_view v_x = gsl_vector_view_array(
                &ptr_x[j * n_parameters], n_parameters);

            gsl_ran_multivariate_gaussian_pdf(&v_xx.vector, &v_x.vector,
                                              SIGMA, &pdf, work);

            sum += ptr_w[j] * pdf;
        }

        sum = log(sum);
        ptr_ww[i] -= sum;
        if (ptr_ww[i] > max_ww)
            max_ww = ptr_ww[i];
    }

    sum = 0.0;
    for (int i = 0; i < n_particles; i++) {
        ptr_ww[i] -= max_ww;
        ptr_ww[i] = exp(ptr_ww[i]);
        sum += ptr_ww[i];
    }

    /* Normalize. */
    for (int i = 0; i < n_particles; i++)
        ptr_ww[i] /= sum;

cleanup:
    gsl_matrix_free(SIGMA);
    gsl_vector_free(work);

    if (error)
        SimInf_abc_error(error);

    UNPROTECT(1);

    return ww;
}

#endif
