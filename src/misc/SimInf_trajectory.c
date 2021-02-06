/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
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

#include <Rdefines.h>
#include <R_ext/Visibility.h>
#include "SimInf.h"
#include "SimInf_openmp.h"
#include "kvec.h"

typedef struct {
    R_xlen_t id;
    R_xlen_t time;
} rowinfo_t;

typedef struct {
    size_t n, m;
    rowinfo_t *a;
} rowinfo_vec;

static void
SimInf_insert_id_time(
    rowinfo_vec *ri,
    SEXP m,
    R_xlen_t m_stride,
    R_xlen_t tlen)
{
    int *m_ir = INTEGER(GET_SLOT(m, Rf_install("i")));
    int *m_jc = INTEGER(GET_SLOT(m, Rf_install("p")));

    for (R_xlen_t t = 0; t < tlen; t++) {
        R_xlen_t id_last = -1;

        for (R_xlen_t j = m_jc[t]; j < m_jc[t + 1]; j++) {
            R_xlen_t id = m_ir[j] / m_stride;

            if (id > id_last) {
                rowinfo_t r = {id, t};
                kv_push(rowinfo_t, *ri, r);
                id_last = id;
            }
        }
    }
}

static void
SimInf_insert_id_time2(
    rowinfo_vec *ri,
    SEXP m1,
    SEXP m2,
    R_xlen_t m1_stride,
    R_xlen_t m2_stride,
    R_xlen_t tlen)
{
    int *m1_ir = INTEGER(GET_SLOT(m1, Rf_install("i")));
    int *m2_ir = INTEGER(GET_SLOT(m2, Rf_install("i")));
    int *m1_jc = INTEGER(GET_SLOT(m1, Rf_install("p")));
    int *m2_jc = INTEGER(GET_SLOT(m2, Rf_install("p")));

    for (R_xlen_t t = 0; t < tlen; t++) {
        R_xlen_t id_last = -1;
        R_xlen_t j1 = m1_jc[t];
        R_xlen_t j2 = m2_jc[t];

        while (j1 < m1_jc[t + 1] || j2 < m2_jc[t + 1]) {
            R_xlen_t id;

            if (j1 < m1_jc[t + 1]) {
                if (j2 < m2_jc[t + 1]) {
                    R_xlen_t id1 = m1_ir[j1] / m1_stride;
                    R_xlen_t id2 = m2_ir[j2] / m2_stride;

                    if (id1 < id2) {
                        id = id1;
                        j1++;
                    } else {
                        id = id2;
                        j2++;
                    }
                } else {
                    id = m1_ir[j1++] / m1_stride;
                }
            } else {
                id = m2_ir[j2++] / m2_stride;
            }

            if (id > id_last) {
                rowinfo_t r = {id, t};
                kv_push(rowinfo_t, *ri, r);
                id_last = id;
            }
        }
    }
}

static void
SimInf_sparse2df_int(
    SEXP dst,
    rowinfo_vec *ri,
    SEXP m,
    int * m_i,
    R_xlen_t m_i_len,
    R_xlen_t m_stride,
    R_xlen_t nrow,
    R_xlen_t tlen,
    R_xlen_t n_id,
    R_xlen_t col)
{
    int *m_ir = INTEGER(GET_SLOT(m, Rf_install("i")));
    int *m_jc = INTEGER(GET_SLOT(m, Rf_install("p")));
    double *m_x = REAL(GET_SLOT(m, Rf_install("x")));

    for (R_xlen_t i = 0; i < m_i_len; i++) {
        SEXP vec = PROTECT(Rf_allocVector(INTSXP, nrow));
        int *p_vec = INTEGER(vec);

        if (ri) {
            size_t k = 0;
            R_xlen_t p_vec_i = 0, j = 0;

            while (k < kv_size(*ri)) {
                R_xlen_t p_time = kv_A(*ri, k).time;

                while (m_jc[p_time] <= j && j < m_jc[p_time + 1]) {
                    /* Check if data for column. */
                    if (m_ir[j] % m_stride == (m_i[i] - 1)) {
                        R_xlen_t m_id = m_ir[j] / m_stride;

                        if (m_id < kv_A(*ri, k).id) {
                            j++; /* Move on. */
                        } else {
                            if (m_id == kv_A(*ri, k).id)
                                p_vec[p_vec_i++] = m_x[j++];
                            else
                                p_vec[p_vec_i++] = NA_INTEGER;

                            if (++k >= kv_size(*ri))
                                break;
                        }
                    } else {
                        j++; /* Move on. */
                    }
                }

                while (k < kv_size(*ri) && kv_A(*ri, k).time <= p_time) {
                    p_vec[p_vec_i++] = NA_INTEGER;
                    k++;
                }
            }
        } else {
            #ifdef _OPENMP
            #  pragma omp parallel for num_threads(SimInf_num_threads())
            #endif
            for (R_xlen_t t = 0; t < tlen; t++) {
                R_xlen_t id = 0;

                for (R_xlen_t j = m_jc[t]; j < m_jc[t + 1]; j++) {
                    if ((m_ir[j] % m_stride) == (m_i[i] - 1)) {
                        R_xlen_t m_id = m_ir[j] / m_stride;

                        for (; id < m_id; id++)
                            p_vec[t * n_id + id] = NA_INTEGER;

                        p_vec[t * n_id + id] = m_x[j];
                        id++;
                    }
                }

                for (; id < n_id; id++)
                    p_vec[t * n_id + id] = NA_INTEGER;
            }
        }

        SET_VECTOR_ELT(dst, col++, vec);
        UNPROTECT(1);
    }
}

static void
SimInf_sparse2df_real(
    SEXP dst,
    rowinfo_vec *ri,
    SEXP m,
    int * m_i,
    R_xlen_t m_i_len,
    R_xlen_t m_stride,
    R_xlen_t nrow,
    R_xlen_t tlen,
    R_xlen_t n_id,
    R_xlen_t col)
{
    int *m_ir = INTEGER(GET_SLOT(m, Rf_install("i")));
    int *m_jc = INTEGER(GET_SLOT(m, Rf_install("p")));
    double *m_x = REAL(GET_SLOT(m, Rf_install("x")));

    for (R_xlen_t i = 0; i < m_i_len; i++) {
        SEXP vec = PROTECT(Rf_allocVector(REALSXP, nrow));
        double *p_vec = REAL(vec);

        if (ri) {
            size_t k = 0;
            R_xlen_t p_vec_i = 0, j = 0;

            while (k < kv_size(*ri)) {
                R_xlen_t p_time = kv_A(*ri, k).time;

                while (m_jc[p_time] <= j && j < m_jc[p_time + 1]) {
                    /* Check if data for column. */
                    if (m_ir[j] % m_stride == (m_i[i] - 1)) {
                        R_xlen_t m_id = m_ir[j] / m_stride;

                        if (m_id < kv_A(*ri, k).id) {
                            j++; /* Move on. */
                        } else {
                            if (m_id == kv_A(*ri, k).id)
                                p_vec[p_vec_i++] = m_x[j++];
                            else
                                p_vec[p_vec_i++] = NA_REAL;

                            if (++k >= kv_size(*ri))
                                break;
                        }
                    } else {
                        j++; /* Move on. */
                    }
                }

                while (k < kv_size(*ri) && kv_A(*ri, k).time <= p_time) {
                    p_vec[p_vec_i++] = NA_REAL;
                    k++;
                }
            }
        } else {
            #ifdef _OPENMP
            #  pragma omp parallel for num_threads(SimInf_num_threads())
            #endif
            for (R_xlen_t t = 0; t < tlen; t++) {
                R_xlen_t id = 0;

                for (R_xlen_t j = m_jc[t]; j < m_jc[t + 1]; j++) {
                    if ((m_ir[j] % m_stride) == (m_i[i] - 1)) {
                        R_xlen_t m_id = m_ir[j] / m_stride;

                        for (; id < m_id; id++)
                            p_vec[t * n_id + id] = NA_REAL;

                        p_vec[t * n_id + id] = m_x[j];
                        id++;
                    }
                }

                for (; id < n_id; id++)
                    p_vec[t * n_id + id] = NA_REAL;
            }
        }

        SET_VECTOR_ELT(dst, col++, vec);
        UNPROTECT(1);
    }
}

static void
SimInf_dense2df_int(
    SEXP dst,
    int *m,
    int * m_i,
    R_xlen_t m_i_len,
    R_xlen_t m_stride,
    R_xlen_t nrow,
    R_xlen_t tlen,
    R_xlen_t id_len,
    R_xlen_t id_n,
    R_xlen_t col,
    int *p_id)
{
    for (R_xlen_t i = 0; i < m_i_len; i++) {
        SEXP vec = PROTECT(Rf_allocVector(INTSXP, nrow));
        int *p_vec = INTEGER(vec);
        int *p_m = m + m_i[i] - 1;

        if (p_id != NULL) {
            /* Note that the identifiers are one-based. */
            #ifdef _OPENMP
            #  pragma omp parallel for num_threads(SimInf_num_threads())
            #endif
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t i = 0; i < id_len; i++) {
                    p_vec[t * id_len + i] =
                        p_m[(t * id_n + p_id[i] - 1) * m_stride];
                }
            }
        } else {
            #ifdef _OPENMP
            #  pragma omp parallel for num_threads(SimInf_num_threads())
            #endif
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t i = 0; i < id_len; i++) {
                    p_vec[t * id_len + i] =
                        p_m[(t * id_n + i) * m_stride];
                }
            }
        }

        SET_VECTOR_ELT(dst, col++, vec);
        UNPROTECT(1);
    }
}

static void
SimInf_dense2df_real(
    SEXP dst,
    double *m,
    int * m_i,
    R_xlen_t m_i_len,
    R_xlen_t m_stride,
    R_xlen_t nrow,
    R_xlen_t tlen,
    R_xlen_t id_len,
    R_xlen_t id_n,
    R_xlen_t col,
    int *p_id)
{
    for (R_xlen_t i = 0; i < m_i_len; i++) {
        SEXP vec = PROTECT(Rf_allocVector(REALSXP, nrow));
        double *p_vec = REAL(vec);
        double *p_m = m + m_i[i] - 1;

        if (p_id != NULL) {
            /* Note that the node identifiers are one-based. */
            #ifdef _OPENMP
            #  pragma omp parallel for num_threads(SimInf_num_threads())
            #endif
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t i = 0; i < id_len; i++) {
                    p_vec[t * id_len + i] =
                        p_m[(t * id_n + p_id[i] - 1) * m_stride];
                }
            }
        } else {
            #ifdef _OPENMP
            #  pragma omp parallel for num_threads(SimInf_num_threads())
            #endif
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t i = 0; i < id_len; i++) {
                    p_vec[t * id_len + i] =
                        p_m[(t * id_n + i) * m_stride];
                }
            }
        }

        SET_VECTOR_ELT(dst, col++, vec);
        UNPROTECT(1);
    }
}

/**
 * Extract data from a simulated trajectory as a data.frame.
 *
 * @param dm data for the discrete state matrix to transform
 *        to a data.frame.
 * @param dm_i index (1-based) to compartments in 'dm' to include
 *        in the data.frame.
 * @param dm_lbl state names of the data in 'dm'.
 * @param cm data for the continuous state matrix to transform
 *        to a data.frame.
 * @param cm_i index (1-based) to compartments in 'cm' to include
 *        in the data.frame.
 * @param cm_lbl state names of the data in 'cm'.
 * @param tspan a vector of increasing time points for the time
 *        in each column in 'dm' and 'cm'.
 * @param id_n number of identifiers in the model.
 * @param id NULL or an integer vector with (1-based) indices of the
 *        identifiers to include in the data.frame.
 * @param id_lbl character vector of length one with the name of the
 *        identifier column.
 * @return A data.frame.
 */
SEXP attribute_hidden
SimInf_trajectory(
    SEXP dm,
    SEXP dm_i,
    SEXP dm_lbl,
    SEXP cm,
    SEXP cm_i,
    SEXP cm_lbl,
    SEXP tspan,
    SEXP id_n,
    SEXP id,
    SEXP id_lbl)
{
    SEXP colnames, result, vec;
    int *p_vec;
    int *p_id = Rf_isNull(id) ? NULL : INTEGER(id);
    R_xlen_t dm_i_len = XLENGTH(dm_i);
    R_xlen_t dm_stride = Rf_isNull(dm_lbl) ? 0 : XLENGTH(dm_lbl);
    int dm_sparse = Rf_isS4(dm) && Rf_inherits(dm, "dgCMatrix") ? 1 : 0;
    R_xlen_t cm_i_len = XLENGTH(cm_i);
    R_xlen_t cm_stride = Rf_isNull(cm_lbl) ? 0 : XLENGTH(cm_lbl);
    int cm_sparse = Rf_isS4(cm) && Rf_inherits(cm, "dgCMatrix") ? 1 : 0;
    R_xlen_t tlen = XLENGTH(tspan);
    R_xlen_t c_id_n = Rf_asInteger(id_n);
    R_xlen_t id_len = Rf_isNull(id) ? c_id_n : XLENGTH(id);
    R_xlen_t nrow = tlen * id_len;
    R_xlen_t ncol = 2 + dm_i_len + cm_i_len; /* The '2' is for the
                                              * 'identifier' and
                                              * 'time' columns. */
    rowinfo_vec *ri = NULL;

    /* Use all available threads in parallel regions. */
    SimInf_set_num_threads(-1);

    /* Create a vector for the column names. */
    PROTECT(colnames = Rf_allocVector(STRSXP, ncol));
    SET_STRING_ELT(colnames, 0, STRING_ELT(id_lbl, 0));
    SET_STRING_ELT(colnames, 1, Rf_mkChar("time"));
    for (R_xlen_t i = 0; i < dm_i_len; i++) {
        R_xlen_t j = INTEGER(dm_i)[i] - 1;
        SET_STRING_ELT(colnames, 2 + i, STRING_ELT(dm_lbl, j));
    }
    for (R_xlen_t i = 0; i < cm_i_len; i++) {
        R_xlen_t j = INTEGER(cm_i)[i] - 1;
        SET_STRING_ELT(colnames, 2 + dm_i_len + i, STRING_ELT(cm_lbl, j));
    }

    /* Determine the number of rows that is required for the
     * data.frame. If either U or V is a dense matrix, then we need a
     * full data.frame with one row per node and time point, else the
     * number of rows depends on unique combinations of identifier and
     * time information in the sparse matrices. */
    if (dm_i_len > 0 && cm_i_len > 0) {
        if (dm_sparse && cm_sparse) {
            ri = calloc(1, sizeof(rowinfo_vec));
            SimInf_insert_id_time2(ri, dm, cm, dm_stride, cm_stride, tlen);
            nrow = kv_size(*ri);
        }
    } else if (dm_i_len > 0 && dm_sparse) {
        ri = calloc(1, sizeof(rowinfo_vec));
        SimInf_insert_id_time(ri, dm, dm_stride, tlen);
        nrow = kv_size(*ri);
    } else if (cm_i_len > 0 && cm_sparse) {
        ri = calloc(1, sizeof(rowinfo_vec));
        SimInf_insert_id_time(ri, cm, cm_stride, tlen);
        nrow = kv_size(*ri);
    }

    /* Create a list for the 'data.frame' and add colnames and a
     * 'data.frame' class attribute. */
    PROTECT(result = Rf_allocVector(VECSXP, ncol));
    Rf_setAttrib(result, R_NamesSymbol, colnames);
    Rf_setAttrib(result, R_ClassSymbol, Rf_mkString("data.frame"));

    /* Add row names to the 'data.frame'. Note that the row names are
     * one-based. */
    PROTECT(vec = Rf_allocVector(INTSXP, nrow));
    p_vec = INTEGER(vec);
    #ifdef _OPENMP
    #  pragma omp parallel for num_threads(SimInf_num_threads())
    #endif
    for (R_xlen_t i = 0; i < nrow; i++) {
        p_vec[i] = i + 1;
    }
    Rf_setAttrib(result, R_RowNamesSymbol, vec);
    UNPROTECT(1);

    /* Add an identifier column to the 'data.frame'. */
    PROTECT(vec = Rf_allocVector(INTSXP, nrow));
    p_vec = INTEGER(vec);
    if (ri) {
        for (size_t i = 0; i < kv_size(*ri); i++)
            p_vec[i] = kv_A(*ri, i).id + 1;
    } else if (p_id != NULL) {
        #ifdef _OPENMP
        #  pragma omp parallel for num_threads(SimInf_num_threads())
        #endif
        for (R_xlen_t t = 0; t < tlen; t++) {
            memcpy(&p_vec[t * id_len], p_id, id_len * sizeof(int));
        }
    } else {
        #ifdef _OPENMP
        #  pragma omp parallel for num_threads(SimInf_num_threads())
        #endif
        for (R_xlen_t t = 0; t < tlen; t++) {
            for (R_xlen_t i = 0; i < id_len; i++)
                p_vec[t * id_len + i] = i + 1;
        }
    }
    SET_VECTOR_ELT(result, 0, vec);
    UNPROTECT(1);

    /* Add a 'time' column to the 'data.frame'. */
    if (Rf_isNull(Rf_getAttrib(tspan, R_NamesSymbol))) {
        double *p_tspan = REAL(tspan);

        PROTECT(vec = Rf_allocVector(INTSXP, nrow));
        p_vec = INTEGER(vec);

        if (ri) {
            for (size_t i = 0; i < kv_size(*ri); i++)
                p_vec[i] = p_tspan[kv_A(*ri, i).time];
        } else {
            #ifdef _OPENMP
            #  pragma omp parallel for num_threads(SimInf_num_threads())
            #endif
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t i = 0; i < id_len; i++)
                    p_vec[t * id_len + i] = p_tspan[t];
            }
        }

        SET_VECTOR_ELT(result, 1, vec);
        UNPROTECT(1);
    } else {
        SEXP lbl_tspan = PROTECT(Rf_getAttrib(tspan, R_NamesSymbol));

        PROTECT(vec = Rf_allocVector(STRSXP, nrow));

        if (ri) {
            for (size_t i = 0; i < kv_size(*ri); i++)
                SET_STRING_ELT(vec, i, STRING_ELT(lbl_tspan, kv_A(*ri, i).time));
        } else {
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t i = 0; i < id_len; i++)
                    SET_STRING_ELT(vec, t * id_len + i, STRING_ELT(lbl_tspan, t));
            }
        }

        SET_VECTOR_ELT(result, 1, vec);
        UNPROTECT(2);
    }

    /* Copy data from the discrete state matrix. */
    if (dm_sparse) {
        SimInf_sparse2df_int(result, ri, dm, INTEGER(dm_i), dm_i_len,
                             dm_stride, nrow, tlen, id_len, 2);
    } else {
        SimInf_dense2df_int(result, INTEGER(dm), INTEGER(dm_i), dm_i_len,
                            dm_stride, nrow, tlen, id_len, c_id_n, 2, p_id);
    }

    /* Copy data from the continuous state matrix. */
    if (cm_sparse) {
        SimInf_sparse2df_real(result, ri, cm, INTEGER(cm_i), cm_i_len,
                              cm_stride, nrow, tlen, id_len, 2 + dm_i_len);
    } else {
        SimInf_dense2df_real(result, REAL(cm), INTEGER(cm_i), cm_i_len, cm_stride,
                             nrow, tlen, id_len, c_id_n, 2 + dm_i_len, p_id);
    }

    if (ri) {
        kv_destroy(*ri);
        free(ri);
    }

    UNPROTECT(2);

    return result;
}
