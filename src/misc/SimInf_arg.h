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

#pragma once

#include <Rinternals.h>

int SimInf_arg_check_dgCMatrix(SEXP arg);
int SimInf_arg_check_integer(SEXP arg);
int SimInf_arg_check_integer_gt_zero(SEXP arg);
int SimInf_arg_check_matrix(SEXP arg);
int SimInf_arg_check_model(SEXP arg);
int SimInf_get_solver(int *out, SEXP solver);
int SimInf_sparse(SEXP m, R_xlen_t i, R_xlen_t j);
