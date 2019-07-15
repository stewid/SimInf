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

/**
 * Extract data from a simulated trajectory as a data.frame.
 *
 * @param model the SimInf_model with data to transform to a
 *        data.frame.
 * @param ui index (1-based) to compartments in 'U' to include
 *        in the data.frame.
 * @param vi index (1-based) to compartments in 'V' to include
 *        in the data.frame.
 * @param nodes NULL or integer vector with (1-based) node
 *        indices of the nodes to include in the data.frame.
 * @return A data.frame.
 */
SEXP SimInf_trajectory(SEXP model, SEXP ui, SEXP vi, SEXP nodes)
{
    int col = 0, *p_int_vec, *p_nodes = NULL;
    double *p_real_vec;
    R_xlen_t ncol, nrow, Nc, Nd, Nn, Nnodes, tlen;
    SEXP result;
    SEXP colnames, S, tspan, vec, U, V, u0, v0;

    /* Use all available threads in parallel regions. */
    SimInf_set_num_threads(-1);

    PROTECT(S = GET_SLOT(model, Rf_install("S")));
    PROTECT(tspan = GET_SLOT(model, Rf_install("tspan")));
    PROTECT(U = GET_SLOT(model, Rf_install("U")));
    PROTECT(V = GET_SLOT(model, Rf_install("V")));
    PROTECT(u0 = GET_SLOT(model, Rf_install("u0")));
    PROTECT(v0 = GET_SLOT(model, Rf_install("v0")));

    Nc = INTEGER(GET_SLOT(S, Rf_install("Dim")))[0];
    Nn = INTEGER(GET_SLOT(u0, R_DimSymbol))[1];
    Nd = INTEGER(GET_SLOT(v0, R_DimSymbol))[0];
    tlen = XLENGTH(tspan);

    if (Rf_isNull(nodes)) {
        Nnodes = Nn;
    } else {
        Nnodes = XLENGTH(nodes);
        p_nodes = INTEGER(nodes);
    }

    /* Determine the dimensions to hold the trajectory data. */
    nrow = tlen * Nnodes;
    ncol = 2 + XLENGTH(ui) + XLENGTH(vi);

    /* Create a list for the 'data.frame'. */
    PROTECT(result = Rf_allocVector(VECSXP, ncol));

    /* Create a vector for the column names. */
    PROTECT(colnames = Rf_allocVector(STRSXP, ncol));
    Rf_setAttrib(result, R_NamesSymbol, colnames);

    /* Add the 'data.frame' class attribute to the list. */
    Rf_setAttrib(result, R_ClassSymbol, Rf_mkString("data.frame"));

    /* Add row names to the 'data.frame'. Note that the row names are
     * one-based. */
    PROTECT(vec = Rf_allocVector(INTSXP, nrow));
    p_int_vec = INTEGER(vec);
    #pragma omp parallel for num_threads(SimInf_num_threads())
    for (R_xlen_t i = 0; i < nrow; i++)
        p_int_vec[i] = i + 1;
    Rf_setAttrib(result, R_RowNamesSymbol, vec);
    UNPROTECT(1);

    /* Add a 'node' identifier column to the 'data.frame'. */
    SET_STRING_ELT(colnames, col, Rf_mkChar("node"));
    PROTECT(vec = Rf_allocVector(INTSXP, nrow));
    p_int_vec = INTEGER(vec);
    if (p_nodes != NULL) {
        #pragma omp parallel for num_threads(SimInf_num_threads())
        for (R_xlen_t t = 0; t < tlen; t++)
            memcpy(&p_int_vec[t * Nnodes], p_nodes, Nnodes * sizeof(int));
    } else {
        #pragma omp parallel for num_threads(SimInf_num_threads())
        for (R_xlen_t t = 0; t < tlen; t++)
            for (R_xlen_t node = 0; node < Nnodes; node++)
                p_int_vec[t * Nnodes + node] = node + 1;
    }
    SET_VECTOR_ELT(result, col++, vec);
    UNPROTECT(1);

    /* Add a 'time' column to the 'data.frame'. */
    SET_STRING_ELT(colnames, col, Rf_mkChar("time"));
    PROTECT(vec = Rf_allocVector(INTSXP, nrow));
    p_int_vec = INTEGER(vec);
    p_real_vec = REAL(tspan);
    #pragma omp parallel for num_threads(SimInf_num_threads())
    for (R_xlen_t t = 0; t < tlen; t++)
        for (R_xlen_t node = 0; node < Nnodes; node++)
            p_int_vec[t * Nnodes + node] = p_real_vec[t];
    SET_VECTOR_ELT(result, col++, vec);
    UNPROTECT(1);

    if (XLENGTH(ui) > 0) {
        SEXP rownames = VECTOR_ELT(GET_SLOT(S, Rf_install("Dimnames")), 0);

        for (R_xlen_t i = 0; i < XLENGTH(ui); i++) {
            R_xlen_t j = INTEGER(ui)[i] - 1;
            int *p_U = INTEGER(U) + j;

            /* Add data for the compartment to the 'data.frame'. */
            SET_STRING_ELT(colnames, col, STRING_ELT(rownames, j));
            PROTECT(vec = Rf_allocVector(INTSXP, nrow));
            p_int_vec = INTEGER(vec);

            if (p_nodes != NULL) {
                /* Note that the node identifiers are one-based. */
                #pragma omp parallel for num_threads(SimInf_num_threads())
                for (R_xlen_t t = 0; t < tlen; t++) {
                    for (R_xlen_t node = 0; node < Nnodes; node++) {
                        p_int_vec[t * Nnodes + node] =
                            p_U[(t * Nn + p_nodes[node] - 1) * Nc];
                    }
                }
            } else {
                #pragma omp parallel for num_threads(SimInf_num_threads())
                for (R_xlen_t t = 0; t < tlen; t++) {
                    for (R_xlen_t node = 0; node < Nnodes; node++) {
                        p_int_vec[t * Nnodes + node] =
                            p_U[(t * Nn + node) * Nc];
                    }
                }
            }

            SET_VECTOR_ELT(result, col++, vec);
            UNPROTECT(1);
        }
    }

    if (XLENGTH(vi) > 0) {
        SEXP rownames = VECTOR_ELT(Rf_getAttrib(v0, R_DimNamesSymbol), 0);

        for (R_xlen_t i = 0; i < XLENGTH(vi); i++) {
            R_xlen_t j = INTEGER(vi)[i] - 1;
            double *p_V = REAL(V) + j;

            /* Add data for the compartment to the 'data.frame'. */
            SET_STRING_ELT(colnames, col, STRING_ELT(rownames, j));
            PROTECT(vec = Rf_allocVector(REALSXP, nrow));
            p_real_vec = REAL(vec);

            if (p_nodes != NULL) {
                /* Note that the node identifiers are one-based. */
                #pragma omp parallel for num_threads(SimInf_num_threads())
                for (R_xlen_t t = 0; t < tlen; t++) {
                    for (R_xlen_t node = 0; node < Nnodes; node++) {
                        p_real_vec[t * Nnodes + node] =
                            p_V[(t * Nn + p_nodes[node] - 1) * Nd];
                    }
                }
            } else {
                #pragma omp parallel for num_threads(SimInf_num_threads())
                for (R_xlen_t t = 0; t < tlen; t++) {
                    for (R_xlen_t node = 0; node < Nnodes; node++) {
                        p_real_vec[t * Nnodes + node] =
                            p_V[(t * Nn + node) * Nd];
                    }
                }
            }

            SET_VECTOR_ELT(result, col++, vec);
            UNPROTECT(1);
        }
    }

    UNPROTECT(8);

    return result;
}
