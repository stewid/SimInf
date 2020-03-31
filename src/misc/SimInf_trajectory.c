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

#include "SimInf.h"
#include "SimInf_openmp.h"
#include "kbtree.h"

#define SIMINF_UNUSED(x) ((void)(x))

typedef struct {
    R_xlen_t id;
    R_xlen_t time;
} rowinfo_t;

static int rowinfo_cmp(rowinfo_t x, rowinfo_t y)
{
    if (x.time < y.time)
        return -1;
    if (x.time > y.time)
        return 1;
    if (x.id < y.id)
        return -1;
    if (x.id > y.id)
        return 1;
    return 0;
}

KBTREE_INIT(rowinfo, rowinfo_t, rowinfo_cmp)

static void SimInf_insert_node_time(
    kbtree_t(rowinfo) *ri,
    SEXP m,
    R_xlen_t m_stride,
    R_xlen_t tlen)
{
    int *m_ir = INTEGER(GET_SLOT(m, Rf_install("i")));
    int *m_jc = INTEGER(GET_SLOT(m, Rf_install("p")));

    for (R_xlen_t t = 0; t < tlen; t++) {
        for (R_xlen_t j = m_jc[t]; j < m_jc[t + 1]; j++) {
            rowinfo_t r = {m_ir[j] / m_stride, t};
            if (!kb_get(rowinfo, ri, r))
                kb_put(rowinfo, ri, r);
        }
    }
}

static void SimInf_sparse2df_int(
    SEXP dst,
    kbtree_t(rowinfo) *ri,
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

        if (ri != NULL) {
            R_xlen_t p_vec_i = 0, j = 0;
            kbitr_t itr;

            kb_itr_first(rowinfo, ri, &itr);
            while (kb_itr_valid(&itr)) {
                rowinfo_t *p = &kb_itr_key(rowinfo_t, &itr);
                R_xlen_t p_time = p->time;

                while (m_jc[p_time] <= j && j < m_jc[p_time + 1]) {
                    /* Check if data for column. */
                    if (m_ir[j] % m_stride == (m_i[i] - 1)) {
                        R_xlen_t m_id = m_ir[j] / m_stride;

                        if (m_id < p->id) {
                            j++; /* Move on. */
                        } else {
                            if (m_id == p->id)
                                p_vec[p_vec_i++] = m_x[j++];
                            else
                                p_vec[p_vec_i++] = NA_INTEGER;

                            kb_itr_next(rowinfo, ri, &itr);
                            if (!kb_itr_valid(&itr))
                                break;
                            p = &kb_itr_key(rowinfo_t, &itr);
                        }
                    } else {
                        j++; /* Move on. */
                    }
                }

                while (kb_itr_valid(&itr) && p->time <= p_time) {
                    p_vec[p_vec_i++] = NA_INTEGER;
                    kb_itr_next(rowinfo, ri, &itr);
                    if (kb_itr_valid(&itr))
                        p = &kb_itr_key(rowinfo_t, &itr);
                }
            }
        } else {
            #pragma omp parallel for num_threads(SimInf_num_threads())
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

static void SimInf_sparse2df_real(
    SEXP dst,
    kbtree_t(rowinfo) *ri,
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

        if (ri != NULL) {
            R_xlen_t p_vec_i = 0, j = 0;
            kbitr_t itr;

            kb_itr_first(rowinfo, ri, &itr);
            while (kb_itr_valid(&itr)) {
                rowinfo_t *p = &kb_itr_key(rowinfo_t, &itr);
                R_xlen_t p_time = p->time;

                while (m_jc[p_time] <= j && j < m_jc[p_time + 1]) {
                    /* Check if data for column. */
                    if (m_ir[j] % m_stride == (m_i[i] - 1)) {
                        R_xlen_t m_id = m_ir[j] / m_stride;

                        if (m_id < p->id) {
                            j++; /* Move on. */
                        } else {
                            if (m_id == p->id)
                                p_vec[p_vec_i++] = m_x[j++];
                            else
                                p_vec[p_vec_i++] = NA_REAL;

                            kb_itr_next(rowinfo, ri, &itr);
                            if (!kb_itr_valid(&itr))
                                break;
                            p = &kb_itr_key(rowinfo_t, &itr);
                        }
                    } else {
                        j++; /* Move on. */
                    }
                }

                while (kb_itr_valid(&itr) && p->time <= p_time) {
                    p_vec[p_vec_i++] = NA_REAL;
                    kb_itr_next(rowinfo, ri, &itr);
                    if (kb_itr_valid(&itr))
                        p = &kb_itr_key(rowinfo_t, &itr);
                }
            }
        } else {
            #pragma omp parallel for num_threads(SimInf_num_threads())
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

static void SimInf_dense2df_int(
    SEXP dst,
    int *m,
    int * m_i,
    R_xlen_t m_i_len,
    R_xlen_t m_stride,
    R_xlen_t nrow,
    R_xlen_t tlen,
    R_xlen_t Nnodes,
    R_xlen_t Nn,
    R_xlen_t col,
    int *p_nodes)
{
    for (R_xlen_t i = 0; i < m_i_len; i++) {
        SEXP vec = PROTECT(Rf_allocVector(INTSXP, nrow));
        int *p_vec = INTEGER(vec);
        int *p_m = m + m_i[i] - 1;

        if (p_nodes != NULL) {
            /* Note that the node identifiers are one-based. */
            #pragma omp parallel for num_threads(SimInf_num_threads())
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t node = 0; node < Nnodes; node++) {
                    p_vec[t * Nnodes + node] =
                        p_m[(t * Nn + p_nodes[node] - 1) * m_stride];
                }
            }
        } else {
            #pragma omp parallel for num_threads(SimInf_num_threads())
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t node = 0; node < Nnodes; node++) {
                    p_vec[t * Nnodes + node] =
                        p_m[(t * Nn + node) * m_stride];
                }
            }
        }

        SET_VECTOR_ELT(dst, col++, vec);
        UNPROTECT(1);
    }
}

static void SimInf_dense2df_real(
    SEXP dst,
    double *m,
    int * m_i,
    R_xlen_t m_i_len,
    R_xlen_t m_stride,
    R_xlen_t nrow,
    R_xlen_t tlen,
    R_xlen_t Nnodes,
    R_xlen_t Nn,
    R_xlen_t col,
    int *p_nodes)
{
    for (R_xlen_t i = 0; i < m_i_len; i++) {
        SEXP vec = PROTECT(Rf_allocVector(REALSXP, nrow));
        double *p_vec = REAL(vec);
        double *p_m = m + m_i[i] - 1;

        if (p_nodes != NULL) {
            /* Note that the node identifiers are one-based. */
            #pragma omp parallel for num_threads(SimInf_num_threads())
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t node = 0; node < Nnodes; node++) {
                    p_vec[t * Nnodes + node] =
                        p_m[(t * Nn + p_nodes[node] - 1) * m_stride];
                }
            }
        } else {
            #pragma omp parallel for num_threads(SimInf_num_threads())
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t node = 0; node < Nnodes; node++) {
                    p_vec[t * Nnodes + node] =
                        p_m[(t * Nn + node) * m_stride];
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
 * @param Nn number of nodes in the SimInf_model.
 * @param nodes NULL or an integer vector with (1-based) node
 *        indices of the nodes to include in the data.frame.
 * @return A data.frame.
 */
SEXP SimInf_trajectory(
    SEXP dm,
    SEXP dm_i,
    SEXP dm_lbl,
    SEXP cm,
    SEXP cm_i,
    SEXP cm_lbl,
    SEXP tspan,
    SEXP Nn,
    SEXP nodes)
{
    SEXP colnames, result, vec;
    int *p_vec;
    int *p_nodes = Rf_isNull(nodes) ? NULL : INTEGER(nodes);
    R_xlen_t dm_i_len = XLENGTH(dm_i);
    R_xlen_t dm_stride = Rf_isNull(dm_lbl) ? 0 : XLENGTH(dm_lbl);
    int dm_sparse = Rf_isS4(dm) && Rf_inherits(dm, "dgCMatrix") ? 1 : 0;
    R_xlen_t cm_i_len = XLENGTH(cm_i);
    R_xlen_t cm_stride = Rf_isNull(cm_lbl) ? 0 : XLENGTH(cm_lbl);
    int cm_sparse = Rf_isS4(cm) && Rf_inherits(cm, "dgCMatrix") ? 1 : 0;
    R_xlen_t tlen = XLENGTH(tspan);
    R_xlen_t c_Nn = Rf_asInteger(Nn);
    R_xlen_t Nnodes = Rf_isNull(nodes) ? c_Nn : XLENGTH(nodes);
    R_xlen_t nrow = tlen * Nnodes;
    R_xlen_t ncol = 2 + dm_i_len + cm_i_len; /* The '2' is for the 'node' and 'time' columns. */
    kbtree_t(rowinfo) *ri = NULL;

    /* Use all available threads in parallel regions. */
    SimInf_set_num_threads(-1);

    /* Create a vector for the column names. */
    PROTECT(colnames = Rf_allocVector(STRSXP, ncol));
    SET_STRING_ELT(colnames, 0, Rf_mkChar("node"));
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
     * number of rows depends on unique combinations of node and time
     * information in the sparse matrices. */
    if (dm_i_len > 0 && cm_i_len > 0) {
        if (dm_sparse && cm_sparse) {
            ri = kb_init(rowinfo, KB_DEFAULT_SIZE);
            SimInf_insert_node_time(ri, dm, dm_stride, tlen);
            SimInf_insert_node_time(ri, cm, cm_stride, tlen);
            nrow = kb_size(ri);
        }
    } else if (dm_i_len > 0 && dm_sparse) {
        ri = kb_init(rowinfo, KB_DEFAULT_SIZE);
        SimInf_insert_node_time(ri, dm, dm_stride, tlen);
        nrow = kb_size(ri);
    } else if (cm_i_len > 0 && cm_sparse) {
        ri = kb_init(rowinfo, KB_DEFAULT_SIZE);
        SimInf_insert_node_time(ri, cm, cm_stride, tlen);
        nrow = kb_size(ri);
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
    #pragma omp parallel for num_threads(SimInf_num_threads())
    for (R_xlen_t ii = 0; ii < nrow; ii++) {
        p_vec[ii] = ii + 1;
    }
    Rf_setAttrib(result, R_RowNamesSymbol, vec);
    UNPROTECT(1);

    /* Add a 'node' identifier column to the 'data.frame'. */
    PROTECT(vec = Rf_allocVector(INTSXP, nrow));
    p_vec = INTEGER(vec);
    if (ri != NULL) {
        kbitr_t itr;

        /* Silence compiler warning unused function. */
        SIMINF_UNUSED(&kb_del_rowinfo);
        SIMINF_UNUSED(&kb_interval_rowinfo);
        SIMINF_UNUSED(&kb_itr_get_rowinfo);

        kb_itr_first(rowinfo, ri, &itr);
        for (R_xlen_t i = 0; kb_itr_valid(&itr); kb_itr_next(rowinfo, ri, &itr)) {
            rowinfo_t *p = &kb_itr_key(rowinfo_t, &itr);
            p_vec[i++] = p->id + 1;
        }
    } else if (p_nodes != NULL) {
        #pragma omp parallel for num_threads(SimInf_num_threads())
        for (R_xlen_t t = 0; t < tlen; t++) {
            memcpy(&p_vec[t * Nnodes], p_nodes, Nnodes * sizeof(int));
        }
    } else {
        #pragma omp parallel for num_threads(SimInf_num_threads())
        for (R_xlen_t t = 0; t < tlen; t++) {
            for (R_xlen_t node = 0; node < Nnodes; node++)
                p_vec[t * Nnodes + node] = node + 1;
        }
    }
    SET_VECTOR_ELT(result, 0, vec);
    UNPROTECT(1);

    /* Add a 'time' column to the 'data.frame'. */
    if (Rf_isNull(Rf_getAttrib(tspan, R_NamesSymbol))) {
        double *p_tspan = REAL(tspan);

        PROTECT(vec = Rf_allocVector(INTSXP, nrow));
        p_vec = INTEGER(vec);

        if (ri != NULL) {
            kbitr_t itr;

            kb_itr_first(rowinfo, ri, &itr);
            for (R_xlen_t i = 0; kb_itr_valid(&itr); kb_itr_next(rowinfo, ri, &itr)) {
                rowinfo_t *p = &kb_itr_key(rowinfo_t, &itr);
                p_vec[i++] = p_tspan[p->time];
            }
        } else {
            #pragma omp parallel for num_threads(SimInf_num_threads())
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t node = 0; node < Nnodes; node++)
                    p_vec[t * Nnodes + node] = p_tspan[t];
            }
        }

        SET_VECTOR_ELT(result, 1, vec);
        UNPROTECT(1);
    } else {
        SEXP lbl_tspan = PROTECT(Rf_getAttrib(tspan, R_NamesSymbol));

        PROTECT(vec = Rf_allocVector(STRSXP, nrow));

        if (ri != NULL) {
            kbitr_t itr;

            kb_itr_first(rowinfo, ri, &itr);
            for (R_xlen_t i = 0; kb_itr_valid(&itr); kb_itr_next(rowinfo, ri, &itr)) {
                rowinfo_t *p = &kb_itr_key(rowinfo_t, &itr);
                SET_STRING_ELT(vec, i++, STRING_ELT(lbl_tspan, p->time));
            }
        } else {
            for (R_xlen_t t = 0; t < tlen; t++) {
                for (R_xlen_t node = 0; node < Nnodes; node++) {
                    SET_STRING_ELT(vec, t * Nnodes + node, STRING_ELT(lbl_tspan, t));
                }
            }
        }

        SET_VECTOR_ELT(result, 1, vec);
        UNPROTECT(2);
    }

    /* Copy data from the discrete state matrix. */
    if (dm_sparse) {
        SimInf_sparse2df_int(result, ri, dm, INTEGER(dm_i), dm_i_len,
                             dm_stride, nrow, tlen, Nnodes, 2);
    } else {
        SimInf_dense2df_int(result, INTEGER(dm), INTEGER(dm_i), dm_i_len,
                            dm_stride, nrow, tlen, Nnodes, c_Nn, 2, p_nodes);
    }

    /* Copy data from the continuous state matrix. */
    if (cm_sparse) {
        SimInf_sparse2df_real(result, ri, cm, INTEGER(cm_i), cm_i_len,
                              cm_stride, nrow, tlen, Nnodes, 2 + dm_i_len);
    } else {
        SimInf_dense2df_real(result, REAL(cm), INTEGER(cm_i), cm_i_len, cm_stride,
                             nrow, tlen, Nnodes, c_Nn, 2 + dm_i_len, p_nodes);
    }

    if (ri != NULL)
        kb_destroy(rowinfo, ri);

    UNPROTECT(2);

    return result;
}
