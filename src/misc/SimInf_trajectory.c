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
 * Get the list element named str, or return NULL.
 */
static SEXP SimInf_get_list_element(SEXP list, const char *str)
{
    R_xlen_t i;
    SEXP elmt = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);

    for (i = 0; i < Rf_xlength(list); i++) {
        if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }
    }

    return elmt;
}

/**
 *  Number of data.frame columns that is required for the trajectory.
 */
static R_xlen_t SimInf_trajectory_ncol(SEXP compartments)
{
    R_xlen_t ncol = 2; /* Identifier and time columns. */
    SEXP elmt;

    elmt = SimInf_get_list_element(compartments, "U");
    if (!Rf_isNull(elmt))
        ncol += XLENGTH(elmt);

    elmt = SimInf_get_list_element(compartments, "V");
    if (!Rf_isNull(elmt))
        ncol += XLENGTH(elmt);

    return ncol;
}

/**
 * Extract data from a simulated trajectory as a data.frame.
 *
 * @param model the SimInf_model with data to transform to a
 *        data.frame.
 * @param compartments a list of character vectors with the
 *        compartments to include in the data.frame.
 * @return A data.frame.
 */
SEXP SimInf_trajectory(SEXP model, SEXP compartments)
{
    int nprotect = 0, col = 0, *p_int_vec;
    double *p_real_vec;
    R_xlen_t ncol, Nc, Nd, Nn, tlen;
    SEXP result;
    SEXP colnames, S, tspan, vec, U, V, u0, v0, elmt;

    /* Use all available threads in parallel regions. */
    SimInf_set_num_threads(-1);

    PROTECT(S = GET_SLOT(model, Rf_install("S")));
    nprotect++;
    PROTECT(tspan = GET_SLOT(model, Rf_install("tspan")));
    nprotect++;
    PROTECT(U = GET_SLOT(model, Rf_install("U")));
    nprotect++;
    PROTECT(V = GET_SLOT(model, Rf_install("V")));
    nprotect++;
    PROTECT(u0 = GET_SLOT(model, Rf_install("u0")));
    nprotect++;
    PROTECT(v0 = GET_SLOT(model, Rf_install("v0")));
    nprotect++;

    Nc = INTEGER(GET_SLOT(S, Rf_install("Dim")))[0];
    Nn = INTEGER(GET_SLOT(u0, R_DimSymbol))[1];
    Nd = INTEGER(GET_SLOT(v0, R_DimSymbol))[0];
    tlen = XLENGTH(tspan);

    /* Create a list for the 'data.frame'. */
    ncol = SimInf_trajectory_ncol(compartments);
    PROTECT(result = Rf_allocVector(VECSXP, ncol));
    nprotect++;

    /* Create a vector for the column names. */
    PROTECT(colnames = Rf_allocVector(STRSXP, ncol));
    nprotect++;
    Rf_setAttrib(result, R_NamesSymbol, colnames);

    /* Add the 'data.frame' class attribute to the list. */
    Rf_setAttrib(result, R_ClassSymbol, Rf_mkString("data.frame"));

    /* Add row names to the 'data.frame'. */
    PROTECT(vec = Rf_allocVector(INTSXP, tlen * Nn));
    nprotect++;
    Rf_setAttrib(result, R_RowNamesSymbol, vec);
    p_int_vec = INTEGER(vec);
    #pragma omp parallel for num_threads(SimInf_num_threads())
    for (R_xlen_t i = 0; i < tlen * Nn; i++)
        p_int_vec[i] = i + 1;

    /* Add a 'node' identifier column to the 'data.frame'. */
    SET_STRING_ELT(colnames, col, Rf_mkChar("node"));
    PROTECT(vec = Rf_allocVector(INTSXP, tlen * Nn));
    nprotect++;
    SET_VECTOR_ELT(result, col++, vec);
    p_int_vec = INTEGER(vec);
    #pragma omp parallel for num_threads(SimInf_num_threads())
    for (R_xlen_t i = 0; i < tlen; i++)
        for (R_xlen_t j = 0; j < Nn; j++)
            p_int_vec[i * Nn + j] = j + 1;

    /* Add a 'time' column to the 'data.frame'. */
    SET_STRING_ELT(colnames, col, Rf_mkChar("time"));
    PROTECT(vec = Rf_allocVector(INTSXP, tlen * Nn));
    nprotect++;
    SET_VECTOR_ELT(result, col++, vec);
    p_int_vec = INTEGER(vec);
    p_real_vec = REAL(tspan);
    #pragma omp parallel for num_threads(SimInf_num_threads())
    for (R_xlen_t i = 0; i < tlen; i++)
        for (R_xlen_t j = 0; j < Nn; j++)
            p_int_vec[i*Nn+j] = p_real_vec[i];

    elmt = SimInf_get_list_element(compartments, "U");
    if (!Rf_isNull(elmt)) {
        SEXP rownames = VECTOR_ELT(GET_SLOT(S, Rf_install("Dimnames")), 0);

        for (R_xlen_t i = 0; i < XLENGTH(elmt); i++) {
            R_xlen_t j;
            int *p_U;

            /* Match compartment in model. */
            for (j = 0; j < Nc; j++) {
                if (strcmp(CHAR(STRING_ELT(elmt, i)),
                           CHAR(STRING_ELT(rownames, j))) == 0) {
                    break;
                }
            }

            if (j >= Nc) {
                Rf_error("Non-existing compartment in model: '%s'.",
                         CHAR(STRING_ELT(elmt, i)));
            }

            /* Add matched compartment column to the 'data.frame'. */
            SET_STRING_ELT(colnames, col, STRING_ELT(elmt, i));
            PROTECT(vec = Rf_allocVector(INTSXP, tlen * Nn));
            nprotect++;
            SET_VECTOR_ELT(result, col++, vec);
            p_int_vec = INTEGER(vec);
            p_U = INTEGER(U) + j;
            #pragma omp parallel for num_threads(SimInf_num_threads())
            for (R_xlen_t k = 0; k < tlen; k++)
                for (R_xlen_t l = 0; l < Nn; l++)
                    p_int_vec[k * Nn + l] = p_U[(k * Nn + l) * Nc];
        }
    }

    elmt = SimInf_get_list_element(compartments, "V");
    if (!Rf_isNull(elmt)) {
        SEXP rownames = VECTOR_ELT(Rf_getAttrib(v0, R_DimNamesSymbol), 0);

        for (R_xlen_t i = 0; i < XLENGTH(elmt); i++) {
            R_xlen_t j;
            double *p_V;

            /* Match compartment in model. */
            for (j = 0; j < Nd; j++) {
                if (strcmp(CHAR(STRING_ELT(elmt, i)),
                           CHAR(STRING_ELT(rownames, j))) == 0) {
                    break;
                }
            }

            if (j >= Nd) {
                Rf_error("Non-existing compartment in model: '%s'.",
                         CHAR(STRING_ELT(elmt, i)));
            }

            /* Add matched compartment column to the 'data.frame'. */
            SET_STRING_ELT(colnames, col, STRING_ELT(elmt, i));
            PROTECT(vec = Rf_allocVector(REALSXP, tlen * Nn));
            nprotect++;
            SET_VECTOR_ELT(result, col++, vec);
            p_real_vec = REAL(vec);
            p_V = REAL(V) + j;
            #pragma omp parallel for num_threads(SimInf_num_threads())
            for (R_xlen_t k = 0; k < tlen; k++)
                for (R_xlen_t l = 0; l < Nn; l++)
                    p_real_vec[k * Nn + l] = p_V[(k * Nn + l) * Nd];
        }
    }

    if (nprotect)
        UNPROTECT(nprotect);

    return result;
}
