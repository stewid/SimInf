/*
 *  siminf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015  Stefan Engblom
 *  Copyright (C) 2015  Stefan Widgren
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
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INCLUDE_SIMINF_ARGS_H
#define INCLUDE_SIMINF_ARGS_H

#include <Rdefines.h>

int siminf_arg_check_dgCMatrix(SEXP arg);
int siminf_arg_check_integer(SEXP arg);
int siminf_arg_check_matrix(SEXP arg);
int siminf_arg_check_real(SEXP arg);
int siminf_arg_check_real_vec(SEXP arg, size_t size);
int siminf_get_seed(unsigned long int *out, SEXP seed);
int siminf_get_threads(int *out, SEXP threads);

#endif
