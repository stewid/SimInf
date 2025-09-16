/*
 * This file is part of SimInf, a framework for stochastic
 * disease spread simulations.
 *
 * Copyright (C) 2015 Pavol Bauer
 * Copyright (C) 2017 -- 2019 Robin Eriksson
 * Copyright (C) 2015 -- 2019 Stefan Engblom
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

#include "SimInf.h"
#include <R_ext/Visibility.h>
#include <Rinternals.h>

/**
 * Check dgCMatrix argument
 *
 * @param arg The arg to check
 * @return 0 if OK, else -1
 */
attribute_hidden int
SimInf_arg_check_dgCMatrix(
    SEXP arg)
{
    if (!Rf_isS4(arg))
        return -1;
    SEXP class_name = Rf_getAttrib(arg, R_ClassSymbol);
    if (0 != strcmp(CHAR(STRING_ELT(class_name, 0)), "dgCMatrix"))
        return -1;
    return 0;
}

/**
 * Check integer argument
 *
 * @param arg The arg to check
 * @return 0 if OK, else -1
 */
attribute_hidden int
SimInf_arg_check_integer(
    SEXP arg)
{
    if (!Rf_isInteger(arg) || Rf_length(arg) != 1
        || NA_INTEGER == INTEGER(arg)[0])
        return -1;
    return 0;
}

/**
 * Check that integer argument is greater than zero.
 *
 * @param arg The arg to check
 * @return 0 if OK, else -1
 */
attribute_hidden int
SimInf_arg_check_integer_gt_zero(
    SEXP arg)
{
    if (SimInf_arg_check_integer(arg))
        return -1;
    if (INTEGER(arg)[0] < 1)
        return -1;
    return 0;
}

/**
 * Check matrix argument
 *
 * @param arg The arg to check
 * @return 0 if OK, else -1
 */
attribute_hidden int
SimInf_arg_check_matrix(
    SEXP arg)
{
    if (!Rf_isMatrix(arg))
        return -1;
    return 0;
}

/**
 * Check model argument
 *
 * @param arg The arg to check
 * @return 0 if OK, else -1
 */
attribute_hidden int
SimInf_arg_check_model(
    SEXP arg)
{
    static const char *valid[] = { "SimInf_model", "" };

    if (!Rf_isS4(arg) || R_check_class_etc(arg, valid) < 0)
        return -1;

    return 0;
}

/**
 * Check if the trajectory data is stored in a sparse matrix.
 *
 * @param m sparse matrix
 * @param i number of rows in the matrix if data is stored
 *        in a sparse matrix.
 * @param j number of columsn in the matrix if data is stored
 *        in a sparse matrix.
 * @return 1 if data is stored in the sparse matrix, else 0.
 */
attribute_hidden int
SimInf_sparse(
    SEXP m,
    ptrdiff_t i,
    ptrdiff_t j)
{
    const int *d = INTEGER(R_do_slot(m, Rf_install("Dim")));
    return d[0] == i && d[1] == j;
}
