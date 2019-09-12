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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "SimInf_openmp.h"

/**
 * Utility function for implementing the Approximate Bayesian
 * Computation Sequential Monte Carlo (ABC-SMC) algorithm of Toni et
 * al. (2009). Samples particles from: first the prior distribution
 * for the first generation, and then from a previous generation of
 * particles.
 *
 * @param distribution character vector with the name of the
 *        distribution for each prior. Each entry must contain one of
 *        'G' (gamma), 'N' (normal) or 'U' (uniform).
 * @param p1 numeric vector with the first hyperparameter for each
 *        prior: G) shape, N) mean, and U) lower bound.
 * @param p2 numeric vector with the second hyperparameter for each
 *        prior: G) rate, N) standard deviation, and U) upper bound.
 * @param n number of proposals to generate.
 * @param x a numeric matrix with a previous generation of particles
 *        or NULL.
 * @param w a numeric vector with weigths for the previous generation
 *        of particles or NULL.
 * @param sigma variance-covariance matrix.
 * @return a numeric matrix with proposals. The matrix also has an
 * attribute 'particle' with an index that indicates which particle it
 * was sampled from.
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
    int error = 0;
    const int k = Rf_length(parameter);
    gsl_rng **rng = NULL;
    gsl_matrix_view v_sigma;
    gsl_matrix *SIGMA = NULL;
    int N = INTEGER(n)[0];
    double *ptr_x = NULL, *ptr_w = NULL, *cumsum_w = NULL;
    double *ptr_p1 = REAL(p1), *ptr_p2 = REAL(p2);
    int len;
    SEXP result, particle, dimnames;
    double *ptr_result;
    int *ptr_particle;
    int Nthread;

    /* Setup result matrix. */
    PROTECT(result = Rf_allocMatrix(REALSXP, k, N));
    PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
    Rf_setAttrib(result, R_DimNamesSymbol, dimnames);
    SET_VECTOR_ELT(dimnames, 0, parameter);
    ptr_result = REAL(result);

    /* Setup vector to record which particle the proposal was sampled
     * from. */
    PROTECT(particle = Rf_allocVector(INTSXP, N));
    Rf_setAttrib(result, Rf_install("particle"), particle);
    ptr_particle = INTEGER(particle);

    /* Specify the number of threads to use. Make sure to not use more
     * threads than the number of particles to sample. */
    Nthread = SimInf_set_num_threads(N);

    /* Setup random number generator. */
    GetRNGstate();
    rng = calloc(Nthread, sizeof(gsl_rng*));
    for (int i = 0; i < Nthread; i++) {
        rng[i] = gsl_rng_alloc(gsl_rng_mt19937);
        if (!rng[i]) {
            error = 1;
            goto cleanup;
        }
        gsl_rng_set(rng[i], unif_rand() * UINT_MAX);
    }
    PutRNGstate();

    if (Rf_isNull(x)) {
        /* First generation: sample from priors. */
        #pragma omp parallel for num_threads(SimInf_num_threads())
        for (int i = 0; i < N; i++) {
            int tr = omp_get_thread_num();

            ptr_particle[i] = NA_INTEGER;
            for (int d = 0; d < k; d++) {
                switch(R_CHAR(STRING_ELT(distribution, d))[0]) {
                case 'G':
                    ptr_result[i * k + d] =
                        gsl_ran_gamma(rng[tr], ptr_p1[d], 1.0 / ptr_p2[d]);
                    break;
                case 'N':
                    ptr_result[i * k + d] =
                        gsl_ran_gaussian(rng[tr], ptr_p2[d]) + ptr_p1[d];
                    break;
                case 'U':
                    ptr_result[i * k + d] =
                        gsl_ran_flat(rng[tr], ptr_p1[d], ptr_p2[d]);
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
    cumsum_w = malloc(len * sizeof(double));
    for (int i = 0; i < len; i++) {
        if (ptr_w[i] < 0) {
            error = 3;
            goto cleanup;
        }
        cumsum_w[i] = ptr_w[i];
        if (i > 0)
            cumsum_w[i] += cumsum_w[i-1];
    }

    #pragma omp parallel for num_threads(SimInf_num_threads())
    for (int i = 0; i < N; i++) {
        int accept;
        gsl_vector_view X;
        gsl_vector_view proposal = gsl_vector_view_array(&ptr_result[i * k], k);
        int tr = omp_get_thread_num();

        do {
            /* Sample a particle from previous generation. Use a
             * binary search to determine the sampled particle based
             * on its weight. */
            int j = 0, j_low = 0, j_high = len - 1;
            double r = gsl_rng_uniform_pos(rng[tr]) * cumsum_w[j_high];
            while (j_high >= j_low) {
                j = (j_low + j_high) / 2;
                if (cumsum_w[j] < r)
                    j_low = j + 1;
                else if (cumsum_w[j] - ptr_w[j] > r)
                    j_high = j - 1;
                else
                    break;
            }
            ptr_particle[i] = j + 1; /* R is one-based. */

            /* Perturbate the particle. */
            X = gsl_vector_view_array(&ptr_x[j * k], k);
            gsl_ran_multivariate_gaussian(rng[tr], &X.vector, SIGMA,
                                          &proposal.vector);

            /* Check that the proposal is valid. */
            accept = 1;
            for (int d = 0; d < k; d++) {
                switch(R_CHAR(STRING_ELT(distribution, d))[0]) {
                case 'G':
                    if (gsl_ran_gamma_pdf(gsl_vector_get(&proposal.vector, d),
                                          ptr_p1[d], 1.0 / ptr_p2[d]) == 0) {
                        accept = 0;
                    }
                    break;
                case 'N':
                    if (gsl_ran_gaussian_pdf(
                            gsl_vector_get(&proposal.vector, d) -
                            ptr_x[j * k + d], ptr_p2[d]) == 0) {
                        accept = 0;
                    }
                    break;
                case 'U':
                    if (gsl_ran_flat_pdf(gsl_vector_get(&proposal.vector, d),
                                         ptr_p1[d], ptr_p2[d]) == 0) {
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
    free(cumsum_w);
    gsl_matrix_free(SIGMA);
    for (int i = 0; i < Nthread; i++)
        gsl_rng_free(rng[i]);
    free(rng);

    if (error) {
        switch (error) {
        case 1:
            Rf_error("Unable to allocate memory buffer.");
            break;
        case 2:
            Rf_error("Unknown distribution.");
            break;
        case 3:
            Rf_error("Negative weight detected.");
            break;
        default:
            Rf_error("Unknown error code: %i.", error);
            break;
        }
    }

    UNPROTECT(3);

    return result;
}
