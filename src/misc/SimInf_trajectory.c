/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015 - 2019  Stefan Widgren
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <Rdefines.h>

#include "SimInf.h"
#include "SimInf_openmp.h"

static void SimInf_dense2df_int(
    SEXP dst,
    int *m,
    int * m_i,
    R_xlen_t m_len,
    R_xlen_t m_stride,
    R_xlen_t nrow,
    R_xlen_t tlen,
    R_xlen_t Nnodes,
    R_xlen_t Nn,
    R_xlen_t col,
    int *p_nodes)
{
    for (R_xlen_t i = 0; i < m_len; i++) {
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
    R_xlen_t m_len,
    R_xlen_t m_stride,
    R_xlen_t nrow,
    R_xlen_t tlen,
    R_xlen_t Nnodes,
    R_xlen_t Nn,
    R_xlen_t col,
    int *p_nodes)
{
    for (R_xlen_t i = 0; i < m_len; i++) {
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
    double *p_tspan = REAL(tspan);
    int *p_nodes = Rf_isNull(nodes) ? NULL : INTEGER(nodes);
    R_xlen_t dm_len = XLENGTH(dm_i);
    R_xlen_t dm_stride = Rf_isNull(dm_lbl) ? 0 : XLENGTH(dm_lbl);
    R_xlen_t cm_len = XLENGTH(cm_i);
    R_xlen_t cm_stride = Rf_isNull(cm_lbl) ? 0 : XLENGTH(cm_lbl);
    R_xlen_t tlen = XLENGTH(tspan);
    R_xlen_t c_Nn = Rf_asInteger(Nn);
    R_xlen_t Nnodes = Rf_isNull(nodes) ? c_Nn : XLENGTH(nodes);
    R_xlen_t nrow = tlen * Nnodes;
    R_xlen_t ncol = 2 + dm_len + cm_len; /* The '2' is for the 'node' and 'time' columns. */

    /* Use all available threads in parallel regions. */
    SimInf_set_num_threads(-1);

    /* Create a vector for the column names. */
    PROTECT(colnames = Rf_allocVector(STRSXP, ncol));
    SET_STRING_ELT(colnames, 0, Rf_mkChar("node"));
    SET_STRING_ELT(colnames, 1, Rf_mkChar("time"));
    for (R_xlen_t i = 0; i < dm_len; i++) {
        R_xlen_t j = INTEGER(dm_i)[i] - 1;
        SET_STRING_ELT(colnames, 2 + i, STRING_ELT(dm_lbl, j));
    }
    for (R_xlen_t i = 0; i < cm_len; i++) {
        R_xlen_t j = INTEGER(cm_i)[i] - 1;
        SET_STRING_ELT(colnames, 2 + dm_len + i, STRING_ELT(cm_lbl, j));
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
    for (R_xlen_t i = 0; i < nrow; i++) {
        p_vec[i] = i + 1;
    }
    Rf_setAttrib(result, R_RowNamesSymbol, vec);
    UNPROTECT(1);

    /* Add a 'node' identifier column to the 'data.frame'. */
    PROTECT(vec = Rf_allocVector(INTSXP, nrow));
    p_vec = INTEGER(vec);
    if (p_nodes != NULL) {
        #pragma omp parallel for num_threads(SimInf_num_threads())
        for (R_xlen_t t = 0; t < tlen; t++) {
            memcpy(&p_vec[t * Nnodes], p_nodes, Nnodes * sizeof(int));
        }
    } else {
        #pragma omp parallel for num_threads(SimInf_num_threads())
        for (R_xlen_t t = 0; t < tlen; t++) {
            for (R_xlen_t node = 0; node < Nnodes; node++) {
                p_vec[t * Nnodes + node] = node + 1;
            }
        }
    }
    SET_VECTOR_ELT(result, 0, vec);
    UNPROTECT(1);

    /* Add a 'time' column to the 'data.frame'. */
    PROTECT(vec = Rf_allocVector(INTSXP, nrow));
    p_vec = INTEGER(vec);
    #pragma omp parallel for num_threads(SimInf_num_threads())
    for (R_xlen_t t = 0; t < tlen; t++) {
        for (R_xlen_t node = 0; node < Nnodes; node++) {
            p_vec[t * Nnodes + node] = p_tspan[t];
        }
    }
    SET_VECTOR_ELT(result, 1, vec);
    UNPROTECT(1);

    /* Copy data from the discrete state matrix. */
    SimInf_dense2df_int(result, INTEGER(dm), INTEGER(dm_i), dm_len, dm_stride,
                        nrow, tlen, Nnodes, c_Nn, 2, p_nodes);

    /* Copy data from the continuous state matrix. */
    SimInf_dense2df_real(result, REAL(cm), INTEGER(cm_i), cm_len, cm_stride,
                         nrow, tlen, Nnodes, c_Nn, 2 + dm_len, p_nodes);

    UNPROTECT(2);

    return result;
}
