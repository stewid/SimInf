/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
 * Copyright (C) 2015 -- 2026 Stefan Widgren
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
#include "SimInf_internal.h"
#include "kvec.h"
#include <R_ext/Visibility.h>
#include <Rdefines.h>
#ifdef _OPENMP
#  include <omp.h>
#endif

/**
 * Record of an identifier and column to extract from a sparse matrix.
 *
 * During trajectory extraction, SimInf iterates over the columns of a
 * sparse result matrix (U_sparse or V_sparse) and records each unique
 * identifier present in each column. Each entry pairs the identifier
 * index with the column index from which it was found.
 *
 * The column index encodes both the time step and the replicate:
 * downstream code can decompose it as:
 *
 *   time      = col % tlen
 *   replicate = col / tlen
 *
 * where tlen is the number of time points in tspan.
 *
 * @var id   Zero-based identifier index (node index) within the model.
 * @var col  Zero-based column index into the sparse matrix, encoding
 *           both time step and replicate.
 */
typedef struct {
    ptrdiff_t id;
    ptrdiff_t col;
} rowinfo_t;

typedef
kvec_t(
    rowinfo_t) rowinfo_vec;

/**
 * Extract identifiers and column indices from a single sparse matrix.
 *
 * Iterates over all columns in the sparse matrix (m) and records each
 * unique identifier found per column into the rowinfo vector (ri).
 * Each entry stores the identifier and the column index. When a
 * subset of identifiers is requested (p_id != NULL), only matching
 * identifiers are included. The column index encodes both the time
 * step and the replicate, which downstream code can decompose using
 * col % tlen and col / tlen.
 *
 * @param ri Output vector to append rowinfo entries to.
 * @param m A dgCMatrix sparse matrix (either U_sparse or V_sparse).
 * @param m_stride Number of matrix rows per identifier (e.g., Nc for U
 *     or Nd for V). Used to compute the identifier from a row index.
 * @param p_id Either NULL to include all identifiers, or a sorted
 *     one-based integer vector of identifiers to include.
 * @param id_len Length of p_id. Ignored if p_id is NULL.
 * @return 0 on success, or SIMINF_ERR_SPARSE_MODEL if m_stride < 1.
 */
static int
SimInf_insert_id_time(
    rowinfo_vec *ri,
    SEXP m,
    const ptrdiff_t m_stride,
    const int *p_id,
    const ptrdiff_t id_len)
{
    const int *m_ir = INTEGER(R_do_slot(m, Rf_install("i")));
    const int *m_jc = INTEGER(R_do_slot(m, Rf_install("p")));
    const ptrdiff_t ncol = INTEGER(R_do_slot(m, Rf_install("Dim")))[1];

    if (m_stride < 1)
        return SIMINF_ERR_SPARSE_MODEL; /* #nocov */

    for (ptrdiff_t col = 0; col < ncol; col++) {
        ptrdiff_t id_last = -1;
        ptrdiff_t i = 0;

        for (ptrdiff_t j = m_jc[col]; j < m_jc[col + 1]; j++) {
            ptrdiff_t id = m_ir[j] / m_stride;

            if (id > id_last) {
                id_last = id;

                if (p_id) {
                    /* Note that the identifiers are one-based. */
                    if (i >= id_len || id < (p_id[i] - 1))
                        continue;
                    i++;
                }

                rowinfo_t r = { id, col };
                kv_push(rowinfo_t, *ri, r);
            }
        }
    }

    return 0;
}

/**
 * Extract identifiers and column indices from two sparse matrices
 * (U_sparse and V_sparse).
 *
 * Merges identifiers from U_sparse and V_sparse, which share the same
 * column structure. For each column, identifiers from both matrices
 * are merged in sorted order, ensuring that each unique identifier
 * appears only once per column in the output. Each entry records the
 * identifier and the column index. When a subset of identifiers is
 * requested (p_id != NULL), only matching identifiers are included.
 *
 * The column index encodes both the time step and the replicate,
 * which downstream code can decompose using col % tlen and col / tlen.
 *
 * @param ri Output vector to append rowinfo entries to.
 * @param u The U_sparse dgCMatrix sparse matrix.
 * @param v The V_sparse dgCMatrix sparse matrix.
 * @param u_stride Number of rows per identifier in U_sparse (e.g.,
 *     Nc).
 * @param v_stride Number of rows per identifier in V_sparse (e.g.,
 *     Nd).
 * @param p_id Either NULL to include all identifiers, or a sorted
 *     one-based integer vector of identifiers to include.
 * @param id_len Length of p_id. Ignored if p_id is NULL.
 * @return 0 on success, or SIMINF_ERR_SPARSE_MODEL if a stride is < 1
 *     or if U_sparse and V_sparse have different numbers of columns.
 */
static int
SimInf_insert_id_time2(
    rowinfo_vec *ri,
    SEXP u,
    SEXP v,
    const ptrdiff_t u_stride,
    const ptrdiff_t v_stride,
    const int *p_id,
    const ptrdiff_t id_len)
{
    const int *u_ir = INTEGER(R_do_slot(u, Rf_install("i")));
    const int *u_jc = INTEGER(R_do_slot(u, Rf_install("p")));
    const ptrdiff_t u_ncol = INTEGER(R_do_slot(u, Rf_install("Dim")))[1];
    const int *v_ir = INTEGER(R_do_slot(v, Rf_install("i")));
    const int *v_jc = INTEGER(R_do_slot(v, Rf_install("p")));
    const ptrdiff_t v_ncol = INTEGER(R_do_slot(v, Rf_install("Dim")))[1];

    if (u_stride < 1 || v_stride < 1 || u_ncol != v_ncol)
        return SIMINF_ERR_SPARSE_MODEL; /* #nocov */

    for (ptrdiff_t col = 0; col < u_ncol; col++) {
        ptrdiff_t id_last = -1;
        ptrdiff_t i = 0;
        ptrdiff_t j1 = u_jc[col];
        ptrdiff_t j2 = v_jc[col];

        while (j1 < u_jc[col + 1] || j2 < v_jc[col + 1]) {
            ptrdiff_t id;

            if (j1 < u_jc[col + 1]) {
                if (j2 < v_jc[col + 1]) {
                    ptrdiff_t id1 = u_ir[j1] / u_stride;
                    ptrdiff_t id2 = v_ir[j2] / v_stride;

                    if (id1 < id2) {
                        id = id1;
                        j1++;
                    } else {
                        id = id2;
                        j2++;
                    }
                } else {
                    id = u_ir[j1++] / u_stride;
                }
            } else {
                id = v_ir[j2++] / v_stride;
            }

            if (id > id_last) {
                id_last = id;

                if (p_id) {
                    /* Note that the identifiers are one-based. */
                    if (i >= id_len || id < (p_id[i] - 1))
                        continue;
                    i++;
                }

                rowinfo_t r = { id, col };
                kv_push(rowinfo_t, *ri, r);
            }
        }
    }

    return 0;
}

/**
 * Create the rowinfo vector for trajectory extraction from sparse matrices.
 *
 * Determines which identifiers and columns need to be extracted from
 * sparse result matrices (U_sparse and/or V_sparse) and populates a
 * rowinfo_vec accordingly. The function dispatches to:
 *
 * - SimInf_insert_id_time2: when both U_sparse and V_sparse are
 *   available, merging identifiers from both matrices.
 * - SimInf_insert_id_time: when only one of U_sparse or V_sparse is
 *   available.
 *
 * If neither matrix is sparse or neither has requested compartments
 * (u_i_len and v_i_len are 0), no rowinfo is created and *out remains
 * NULL.
 *
 * @param out Output pointer to a newly allocated rowinfo_vec. Set to
 *     NULL if no sparse extraction is needed.
 * @param u The discrete state result (U or U_sparse).
 * @param v The continuous state result (V or V_sparse).
 * @param u_i_len Number of compartments requested from U.
 * @param v_i_len Number of compartments requested from V.
 * @param u_sparse Non-zero if u is a dgCMatrix sparse matrix.
 * @param v_sparse Non-zero if v is a dgCMatrix sparse matrix.
 * @param u_stride Number of rows per identifier in U_sparse.
 * @param v_stride Number of rows per identifier in V_sparse.
 * @param p_id Either NULL to include all identifiers, or a sorted
 *     one-based integer vector of identifiers to include.
 * @param id_len Length of p_id. Ignored if p_id is NULL.
 * @return 0 on success, or SIMINF_ERR_ALLOC_MEMORY_BUFFER if
 *     allocation fails.
 */
static int
SimInf_create_rowinfo(
    rowinfo_vec **out,
    SEXP u,
    SEXP v,
    const ptrdiff_t u_i_len,
    const ptrdiff_t v_i_len,
    const int u_sparse,
    const int v_sparse,
    const ptrdiff_t u_stride,
    const ptrdiff_t v_stride,
    const int *p_id,
    const ptrdiff_t id_len)
{
    if (u_i_len > 0 && v_i_len > 0) {
        if (u_sparse && v_sparse) {
            *out = calloc(1, sizeof(rowinfo_vec));
            if (!*out)
                return SIMINF_ERR_ALLOC_MEMORY_BUFFER;      /* #nocov */

            return SimInf_insert_id_time2(*out, u, v, u_stride,
                                          v_stride, p_id, id_len);
        }
    } else if (u_i_len > 0 && u_sparse) {
        *out = calloc(1, sizeof(rowinfo_vec));
        if (!*out)
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;          /* #nocov */

        return SimInf_insert_id_time(*out, u, u_stride, p_id, id_len);
    } else if (v_i_len > 0 && v_sparse) {
        *out = calloc(1, sizeof(rowinfo_vec));
        if (!*out)
            return SIMINF_ERR_ALLOC_MEMORY_BUFFER;          /* #nocov */

        return SimInf_insert_id_time(*out, v, v_stride, p_id, id_len);
    }

    return 0;
}

/**
 * Calculate the number of rows required for the trajectory output.
 *
 * Returns the appropriate row count for constructing the result
 * data.frame in trajectory extraction. The calculation differs
 * depending on whether sparse or dense result matrices are used:
 *
 * - Sparse: The number of unique identifier-column combinations found
 *   in the sparse matrix is stored in the rowinfo vector (ri). Return
 *   the size of ri directly.
 *
 * - Dense: The output contains one row per combination of time step,
 *   identifier, and replicate. Return tlen × id_len × replicates.
 *
 * @param ri Pointer to the rowinfo vector populated by
 *     SimInf_create_rowinfo, or NULL if no sparse extraction is
 *     needed.
 * @param tlen Number of time points in tspan.
 * @param id_len Number of identifiers (nodes) to include in the
 *     output.
 * @param replicates Number of model replicates.
 * @return Number of rows for the trajectory output data.frame.
 */
static ptrdiff_t
SimInf_number_of_rows(
    const rowinfo_vec *ri,
    const ptrdiff_t tlen,
    const ptrdiff_t id_len,
    const ptrdiff_t replicates)
{
    if (ri)
        return kv_size(*ri);
    return tlen * id_len * replicates;
}

/**
 * Extract discrete state data from a sparse matrix into a data.frame.
 *
 * Populates one column per requested compartment in the output
 * data.frame destination (dst), starting at dst_col. The function
 * operates in two modes:
 *
 * - With rowinfo (ri != NULL): Iterates over the rowinfo entries,
 *   matching sparse data values to the corresponding identifier and
 *   column. Entries without data are filled with NA_INTEGER. This
 *   mode is used when trajectory extraction is filtered to specific
 *   identifiers.
 *
 * - Without rowinfo (ri == NULL): Iterates over all columns in the
 *   sparse matrix and fills data for all identifiers. Missing entries
 *   are filled with NA_INTEGER. This mode extracts the full matrix
 *   content.
 *
 * @param dst The output list (data.frame columns) to populate.
 * @param ri Rowinfo vector, or NULL for full extraction.
 * @param u The U_sparse dgCMatrix sparse matrix.
 * @param u_i One-based indices of compartments to extract from U.
 * @param u_i_len Number of compartments to extract.
 * @param u_stride Number of rows per identifier in U_sparse (Nc).
 * @param nrow Number of rows in the output columns.
 * @param n_id Number of identifiers (nodes) in the output.
 * @param dst_col Starting column offset in dst for the extracted data.
 */
static void
SimInf_u_sparse2df(
    SEXP dst,
    rowinfo_vec *ri,
    SEXP u,
    const int *u_i,
    const ptrdiff_t u_i_len,
    const ptrdiff_t u_stride,
    const ptrdiff_t nrow,
    const ptrdiff_t n_id,
    const ptrdiff_t dst_col)
{
    const int *u_ir = INTEGER(R_do_slot(u, Rf_install("i")));
    const int *u_jc = INTEGER(R_do_slot(u, Rf_install("p")));
    const double *u_x = REAL(R_do_slot(u, Rf_install("x")));
    const ptrdiff_t ncol = INTEGER(R_do_slot(u, Rf_install("Dim")))[1];

    for (ptrdiff_t i = 0; i < u_i_len; i++) {
        SEXP vec;
        SET_VECTOR_ELT(dst, dst_col + i, vec = Rf_allocVector(INTSXP, nrow));
        int *p_vec = INTEGER(vec);

        if (ri) {
            size_t k = 0;
            ptrdiff_t p_vec_i = 0, j = 0;

            while (k < kv_size(*ri)) {
                const ptrdiff_t col = kv_A(*ri, k).col;

                while (u_jc[col] <= j && j < u_jc[col + 1]) {
                    /* Check if data for column. */
                    if (u_ir[j] % u_stride == (u_i[i] - 1)) {
                        ptrdiff_t u_id = u_ir[j] / u_stride;

                        if (u_id < kv_A(*ri, k).id) {
                            j++;        /* Move on. */
                        } else {
                            if (u_id == kv_A(*ri, k).id)
                                p_vec[p_vec_i++] = (int) u_x[j++];
                            else
                                p_vec[p_vec_i++] = NA_INTEGER;

                            if (++k >= kv_size(*ri))
                                break;
                        }
                    } else {
                        j++;    /* Move on. */
                    }
                }

                while (k < kv_size(*ri) && kv_A(*ri, k).col <= col) {
                    p_vec[p_vec_i++] = NA_INTEGER;
                    k++;
                }
            }
        } else {
            for (ptrdiff_t col = 0; col < ncol; col++) {
                ptrdiff_t id = 0;

                for (ptrdiff_t j = u_jc[col]; j < u_jc[col + 1]; j++) {
                    if ((u_ir[j] % u_stride) == (u_i[i] - 1)) {
                        ptrdiff_t u_id = u_ir[j] / u_stride;

                        for (; id < u_id; id++)
                            p_vec[col * n_id + id] = NA_INTEGER;

                        p_vec[col * n_id + id] = (int) u_x[j];
                        id++;
                    }
                }

                for (; id < n_id; id++)
                    p_vec[col * n_id + id] = NA_INTEGER;
            }
        }
    }
}

/**
 * Extract continuous state data from a sparse matrix into a data.frame.
 *
 * Populates one column per requested continuous state variable in the
 * output data.frame destination (dst), starting at dst_col. The
 * function operates in two modes:
 *
 * - With rowinfo (ri != NULL): Iterates over the rowinfo entries,
 *   matching sparse data values to the corresponding identifier and
 *   column. Entries without data are filled with NA_REAL. This mode
 *   is used when trajectory extraction is filtered to specific
 *   identifiers.
 *
 * - Without rowinfo (ri == NULL): Iterates over all columns in the
 *   sparse matrix and fills data for all identifiers. Missing entries
 *   are filled with NA_REAL. This mode extracts the full matrix
 *   content.
 *
 * @param dst The output list (data.frame columns) to populate.
 * @param ri Rowinfo vector, or NULL for full extraction.
 * @param v The V_sparse dgCMatrix sparse matrix.
 * @param v_i One-based indices of continuous variables to extract from V.
 * @param v_i_len Number of continuous variables to extract.
 * @param v_stride Number of rows per identifier in V_sparse (Nd).
 * @param nrow Number of rows in the output columns.
 * @param n_id Number of identifiers (nodes) in the output.
 * @param dst_col Starting column offset in dst for the extracted data.
 */
static void
SimInf_v_sparse2df(
    SEXP dst,
    rowinfo_vec *ri,
    SEXP v,
    const int *v_i,
    const ptrdiff_t v_i_len,
    const ptrdiff_t v_stride,
    const ptrdiff_t nrow,
    const ptrdiff_t n_id,
    const ptrdiff_t dst_col)
{
    const int *v_ir = INTEGER(R_do_slot(v, Rf_install("i")));
    const int *v_jc = INTEGER(R_do_slot(v, Rf_install("p")));
    const double *v_x = REAL(R_do_slot(v, Rf_install("x")));
    const ptrdiff_t ncol = INTEGER(R_do_slot(v, Rf_install("Dim")))[1];

    for (ptrdiff_t i = 0; i < v_i_len; i++) {
        SEXP vec;
        SET_VECTOR_ELT(dst, dst_col + i, vec = Rf_allocVector(REALSXP, nrow));
        double *p_vec = REAL(vec);

        if (ri) {
            size_t k = 0;
            ptrdiff_t p_vec_i = 0, j = 0;

            while (k < kv_size(*ri)) {
                const ptrdiff_t col = kv_A(*ri, k).col;

                while (v_jc[col] <= j && j < v_jc[col + 1]) {
                    /* Check if data for column. */
                    if (v_ir[j] % v_stride == (v_i[i] - 1)) {
                        ptrdiff_t v_id = v_ir[j] / v_stride;

                        if (v_id < kv_A(*ri, k).id) {
                            j++;        /* Move on. */
                        } else {
                            if (v_id == kv_A(*ri, k).id)
                                p_vec[p_vec_i++] = v_x[j++];
                            else
                                p_vec[p_vec_i++] = NA_REAL;

                            if (++k >= kv_size(*ri))
                                break;
                        }
                    } else {
                        j++;    /* Move on. */
                    }
                }

                while (k < kv_size(*ri) && kv_A(*ri, k).col <= col) {
                    p_vec[p_vec_i++] = NA_REAL;
                    k++;
                }
            }
        } else {
            for (ptrdiff_t col = 0; col < ncol; col++) {
                ptrdiff_t id = 0;

                for (ptrdiff_t j = v_jc[col]; j < v_jc[col + 1]; j++) {
                    if ((v_ir[j] % v_stride) == (v_i[i] - 1)) {
                        ptrdiff_t v_id = v_ir[j] / v_stride;

                        for (; id < v_id; id++)
                            p_vec[col * n_id + id] = NA_REAL;

                        p_vec[col * n_id + id] = v_x[j];
                        id++;
                    }
                }

                for (; id < n_id; id++)
                    p_vec[col * n_id + id] = NA_REAL;
            }
        }
    }
}

/**
 * Extract discrete state data from a dense matrix into a data.frame.
 *
 * Populates one column per requested compartment in the output
 * data.frame destination (dst), starting at dst_col. The function
 * operates in two modes:
 *
 * - With identifier selection (p_id != NULL): Extracts data only for
 *   the requested identifiers (nodes), using their one-based indices
 *   from p_id. The output is ordered by time, then replicate, then
 *   selected identifier.
 *
 * - Without identifier selection (p_id == NULL): Extracts data for all
 *   identifiers in order. The output is ordered by time, then replicate,
 *   then identifier.
 *
 * @param dst The output list (data.frame columns) to populate.
 * @param u Pointer to the dense U matrix data (integer).
 * @param u_i One-based indices of compartments to extract from U.
 * @param u_i_len Number of compartments to extract.
 * @param u_stride Number of rows per identifier in U (Nc).
 * @param nrow Number of rows in the output columns.
 * @param tlen Number of time points in tspan.
 * @param id_len Number of identifiers (nodes) in the output.
 * @param id_n Total number of identifiers (nodes) in the model.
 * @param dst_col Starting column offset in dst for the extracted data.
 * @param p_id Either NULL to extract all identifiers, or a sorted
 *     one-based integer vector of identifiers to include.
 * @param replicates Number of model replicates.
 */
static void
SimInf_u_dense2df(
    SEXP dst,
    const int *u,
    const int *u_i,
    const ptrdiff_t u_i_len,
    const ptrdiff_t u_stride,
    const ptrdiff_t nrow,
    const ptrdiff_t tlen,
    const ptrdiff_t id_len,
    const ptrdiff_t id_n,
    const ptrdiff_t dst_col,
    const int *p_id,
    const ptrdiff_t replicates)
{
    for (ptrdiff_t i = 0; i < u_i_len; i++) {
        const int *p_u = u + u_i[i] - 1;

        SEXP vec;
        SET_VECTOR_ELT(dst, dst_col + i, vec = Rf_allocVector(INTSXP, nrow));
        int *p_vec = INTEGER(vec);

        if (p_id != NULL) {
            /* Note that the identifiers are one-based. */
#ifdef _OPENMP
#  pragma omp parallel for num_threads(SimInf_num_threads())
#endif
            for (ptrdiff_t t = 0; t < tlen; t++) {
                const ptrdiff_t j1 = t * id_len;
                const ptrdiff_t j2 = t * id_n;
                for (ptrdiff_t r = 0; r < replicates; r++) {
                    const ptrdiff_t k1 = r * tlen * id_len;
                    const ptrdiff_t k2 = r * tlen * id_n;
                    for (ptrdiff_t l = 0; l < id_len; l++) {
                        p_vec[j1 + k1 + l] =
                            p_u[(j2 + k2 + p_id[l] - 1) * u_stride];
                    }
                }
            }
        } else {
#ifdef _OPENMP
#  pragma omp parallel for num_threads(SimInf_num_threads())
#endif
            for (ptrdiff_t t = 0; t < tlen; t++) {
                const ptrdiff_t j = t * id_len;
                for (ptrdiff_t r = 0; r < replicates; r++) {
                    const ptrdiff_t k = r * tlen * id_len;
                    for (ptrdiff_t l = 0; l < id_len; l++) {
                        p_vec[j + k + l] = p_u[(j + k + l) * u_stride];
                    }
                }
            }
        }
    }
}

/**
 * Extract continuous state data from a dense matrix into a data.frame.
 *
 * Populates one column per requested continuous state variable in the
 * output data.frame destination (dst), starting at dst_col. The
 * function operates in two modes:
 *
 * - With identifier selection (p_id != NULL): Extracts data only for
 *   the requested identifiers (nodes), using their one-based indices
 *   from p_id. The output is ordered by time, then replicate, then
 *   selected identifier.
 *
 * - Without identifier selection (p_id == NULL): Extracts data for all
 *   identifiers in order. The output is ordered by time, then replicate,
 *   then identifier.
 *
 * @param dst The output list (data.frame columns) to populate.
 * @param v Pointer to the dense V matrix data (double).
 * @param v_i One-based indices of continuous variables to extract from V.
 * @param v_i_len Number of continuous variables to extract.
 * @param v_stride Number of rows per identifier in V (Nd).
 * @param nrow Number of rows in the output columns.
 * @param tlen Number of time points in tspan.
 * @param id_len Number of identifiers (nodes) in the output.
 * @param id_n Total number of identifiers (nodes) in the model.
 * @param dst_col Starting column offset in dst for the extracted data.
 * @param p_id Either NULL to extract all identifiers, or a sorted
 *     one-based integer vector of identifiers to include.
 * @param replicates Number of model replicates.
 */
static void
SimInf_v_dense2df(
    SEXP dst,
    const double *v,
    const int *v_i,
    const ptrdiff_t v_i_len,
    const ptrdiff_t v_stride,
    const ptrdiff_t nrow,
    const ptrdiff_t tlen,
    const ptrdiff_t id_len,
    const ptrdiff_t id_n,
    const ptrdiff_t dst_col,
    const int *p_id,
    const ptrdiff_t replicates)
{
    for (ptrdiff_t i = 0; i < v_i_len; i++) {
        const double *p_v = v + v_i[i] - 1;

        SEXP vec;
        SET_VECTOR_ELT(dst, dst_col + i, vec = Rf_allocVector(REALSXP, nrow));
        double *p_vec = REAL(vec);

        if (p_id != NULL) {
            /* Note that the node identifiers are one-based. */
#ifdef _OPENMP
#  pragma omp parallel for num_threads(SimInf_num_threads())
#endif
            for (ptrdiff_t t = 0; t < tlen; t++) {
                const ptrdiff_t j1 = t * id_len;
                const ptrdiff_t j2 = t * id_n;
                for (ptrdiff_t r = 0; r < replicates; r++) {
                    const ptrdiff_t k1 = r * tlen * id_len;
                    const ptrdiff_t k2 = r * tlen * id_n;
                    for (ptrdiff_t l = 0; l < id_len; l++) {
                        p_vec[j1 + k1 + l] =
                            p_v[(j2 + k2 + p_id[l] - 1) * v_stride];
                    }
                }
            }
        } else {
#ifdef _OPENMP
#  pragma omp parallel for num_threads(SimInf_num_threads())
#endif
            for (ptrdiff_t t = 0; t < tlen; t++) {
                const ptrdiff_t j = t * id_len;
                for (ptrdiff_t r = 0; r < replicates; r++) {
                    const ptrdiff_t k = r * tlen * id_len;
                    for (ptrdiff_t l = 0; l < id_len; l++) {
                        p_vec[j + k + l] = p_v[(j + k + l) * v_stride];
                    }
                }
            }
        }
    }
}

/**
 * Extract data from a simulated trajectory as a data.frame.
 *
 * @param u data for the discrete state matrix to transform to a
 *        data.frame.
 * @param u_i index (1-based) to compartments in 'u' to include in the
 *        data.frame.
 * @param u_lbl state names of the data in 'u'.
 * @param v data for the continuous state matrix to transform to a
 *        data.frame.
 * @param v_i index (1-based) to compartments in 'v' to include in the
 *        data.frame.
 * @param v_lbl state names of the data in 'v'.
 * @param tspan a vector of increasing time points for the time
 *        in each column in 'u' and 'v'.
 * @param id_n number of identifiers in the model.
 * @param id NULL or an integer vector with (1-based) indices of the
 *        identifiers to include in the data.frame.
 * @param id_lbl character vector of length one with the name of the
 *        identifier column.
 * @param n_replicates the number of replicates in the model. This is
 *        only used when u and/or v are dense matrices.
 * @return A data.frame.
 */
attribute_hidden SEXP
SimInf_trajectory(
    SEXP u,
    SEXP u_i,
    SEXP u_lbl,
    SEXP v,
    SEXP v_i,
    SEXP v_lbl,
    SEXP tspan,
    SEXP id_n,
    SEXP id,
    SEXP id_lbl,
    SEXP n_replicates)
{
    int err = 0;
    int nprotect = 0;
    const int *p_id = Rf_isNull(id) ? NULL : INTEGER(id);
    const R_xlen_t u_i_len = XLENGTH(u_i);
    const R_xlen_t u_stride = Rf_isNull(u_lbl) ? 0 : XLENGTH(u_lbl);
    const int u_sparse = Rf_isS4(u) && Rf_inherits(u, "dgCMatrix") ? 1 : 0;
    const R_xlen_t v_i_len = XLENGTH(v_i);
    const R_xlen_t v_stride = Rf_isNull(v_lbl) ? 0 : XLENGTH(v_lbl);
    const int v_sparse = Rf_isS4(v) && Rf_inherits(v, "dgCMatrix") ? 1 : 0;
    const R_xlen_t tlen = XLENGTH(tspan);
    const R_xlen_t c_id_n = Rf_asInteger(id_n);
    const R_xlen_t id_len = Rf_isNull(id) ? c_id_n : XLENGTH(id);
    const R_xlen_t replicates = Rf_asInteger(n_replicates);

    /* Determine the number of columns in the resulting
     * data.frame. The '2' is for the identifier' and 'time'
     * columns. The '(replicates > 1)' is to add a column for the
     * replicate of an identifer and time. */
    const R_xlen_t ncol = 2 + (replicates > 1) + u_i_len + v_i_len;

    /* Use all available threads in parallel regions. */
    SimInf_set_num_threads(-1);

    /* Create a vector for the column names. */
    SEXP colnames = PROTECT(Rf_allocVector(STRSXP, ncol));
    nprotect++;
    R_xlen_t col = 0;
    SET_STRING_ELT(colnames, col++, STRING_ELT(id_lbl, 0));
    SET_STRING_ELT(colnames, col++, Rf_mkChar("time"));
    if (replicates > 1)
        SET_STRING_ELT(colnames, col++, Rf_mkChar("replicate"));
    for (ptrdiff_t i = 0; i < u_i_len; i++) {
        const R_xlen_t j = INTEGER(u_i)[i] - 1;
        SET_STRING_ELT(colnames, col++, STRING_ELT(u_lbl, j));
    }
    for (ptrdiff_t i = 0; i < v_i_len; i++) {
        const R_xlen_t j = INTEGER(v_i)[i] - 1;
        SET_STRING_ELT(colnames, col++, STRING_ELT(v_lbl, j));
    }

    /* Determine the number of rows that is required for the
     * data.frame. If either U or V is a dense matrix, then we need a
     * full data.frame with one row per node and time point, else the
     * number of rows depends on unique combinations of identifier and
     * time information in the sparse matrices. */
    rowinfo_vec *ri = NULL;
    err = SimInf_create_rowinfo(&ri, u, v, u_i_len, v_i_len, u_sparse,
                                v_sparse, u_stride, v_stride, p_id, id_len);
    if (err)
        goto cleanup;           /* #nocov */
    const R_xlen_t nrow = SimInf_number_of_rows(ri, tlen, id_len, replicates);

    /* Create a list for the 'data.frame' and add colnames and a
     * 'data.frame' class attribute. */
    SEXP result = PROTECT(Rf_allocVector(VECSXP, ncol));
    nprotect++;
    Rf_setAttrib(result, R_NamesSymbol, colnames);
    Rf_setAttrib(result, R_ClassSymbol, Rf_mkString("data.frame"));

    /* Add row names to the 'data.frame'. Note that the row names are
     * one-based. */
    SEXP vec = PROTECT(Rf_allocVector(INTSXP, nrow));
    nprotect++;
    int *p_vec = INTEGER(vec);
#ifdef _OPENMP
#  pragma omp parallel for num_threads(SimInf_num_threads())
#endif
    for (ptrdiff_t i = 0; i < nrow; i++) {
        p_vec[i] = (int) (i + 1);
    }
    Rf_setAttrib(result, R_RowNamesSymbol, vec);

    /* Add an identifier column to the 'data.frame'. */
    col = 0;
    SET_VECTOR_ELT(result, col++, vec = Rf_allocVector(INTSXP, nrow));
    p_vec = INTEGER(vec);
    if (ri) {
        for (size_t i = 0; i < kv_size(*ri); i++)
            p_vec[i] = kv_A(*ri, i).id + 1;
    } else {
        for (ptrdiff_t t = 0; t < tlen; t++) {
            const ptrdiff_t j = t * id_len;
            for (ptrdiff_t r = 0; r < replicates; r++) {
                const ptrdiff_t k = r * tlen * id_len;
                if (p_id) {
                    memcpy(&p_vec[j + k], p_id, id_len * sizeof(int));
                } else {
                    for (ptrdiff_t l = 0; l < id_len; l++)
                        p_vec[j + k + l] = (int) (l + 1);
                }
            }
        }
    }

    /* Add a 'time' column to the 'data.frame'. */
    if (Rf_isNull(Rf_getAttrib(tspan, R_NamesSymbol))) {
        const double *p_tspan = REAL(tspan);

        SET_VECTOR_ELT(result, col++, vec = Rf_allocVector(INTSXP, nrow));
        p_vec = INTEGER(vec);
        if (ri) {
            for (size_t i = 0; i < kv_size(*ri); i++)
                p_vec[i] = (int) p_tspan[kv_A(*ri, i).col % tlen];
        } else {
            for (ptrdiff_t t = 0; t < tlen; t++) {
                const ptrdiff_t j = t * id_len;
                for (ptrdiff_t r = 0; r < replicates; r++) {
                    const ptrdiff_t k = r * tlen * id_len;
                    for (ptrdiff_t l = 0; l < id_len; l++)
                        p_vec[j + k + l] = (int) p_tspan[t];
                }
            }
        }
    } else {
        SEXP lbl_tspan = PROTECT(Rf_getAttrib(tspan, R_NamesSymbol));
        nprotect++;

        SET_VECTOR_ELT(result, col++, vec = Rf_allocVector(STRSXP, nrow));
        if (ri) {
            for (ptrdiff_t i = 0; i < kv_size(*ri); i++) {
                SET_STRING_ELT(vec, i,
                               STRING_ELT(lbl_tspan, kv_A(*ri, i).col % tlen));
            }
        } else {
            for (ptrdiff_t t = 0; t < tlen; t++) {
                for (ptrdiff_t i = 0; i < id_len; i++) {
                    SET_STRING_ELT(vec, t * id_len + i,
                                   STRING_ELT(lbl_tspan, t));
                }
            }
        }
    }

    if (replicates > 1) {
        SET_VECTOR_ELT(result, col++, vec = Rf_allocVector(INTSXP, nrow));
        p_vec = INTEGER(vec);

        if (ri) {
            for (size_t i = 0; i < kv_size(*ri); i++)
                p_vec[i] = (int) kv_A(*ri, i).col / tlen + 1;
        } else {
            for (ptrdiff_t r = 0; r < replicates; r++) {
                const ptrdiff_t n = tlen * id_len;
                const ptrdiff_t j = r * n;
                for (ptrdiff_t i = 0; i < n; i++)
                    p_vec[j + i] = (int) (r + 1);
            }
        }
    }

    /* Copy data from the discrete state matrix. */
    if (u_sparse) {
        SimInf_u_sparse2df(result, ri, u, INTEGER(u_i), u_i_len,
                           u_stride, nrow, id_len, col);
    } else {
        SimInf_u_dense2df(result, INTEGER(u), INTEGER(u_i), u_i_len,
                          u_stride, nrow, tlen, id_len, c_id_n, col, p_id,
                          replicates);
    }

    /* Copy data from the continuous state matrix. */
    col += u_i_len;
    if (v_sparse) {
        SimInf_v_sparse2df(result, ri, v, INTEGER(v_i), v_i_len,
                           v_stride, nrow, id_len, col);
    } else {
        SimInf_v_dense2df(result, REAL(v), INTEGER(v_i), v_i_len,
                          v_stride, nrow, tlen, id_len, c_id_n, col, p_id,
                          replicates);
    }

cleanup:
    if (ri) {
        kv_destroy(*ri);
        free(ri);
    }

    if (nprotect)
        UNPROTECT(nprotect);

    if (err)
        Rf_error("Unable to allocate memory buffer.");  /* #nocov */

    return result;
}
