## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2025 Stefan Widgren
##
## SimInf is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SimInf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(SimInf)

## Specify the number of threads to use.
set_num_threads(1)

## For debugging
sessionInfo()

## Define a tolerance
tol <- 1e-8

stopifnot(identical(
    .Call(SimInf:::SimInf_lambertW0, c(NA_real_, Inf, -Inf, NaN)),
    c(NA_real_, Inf, NaN, NaN)))

stopifnot(all(abs(.Call(SimInf:::SimInf_lambertW0, 1) - 0.5671433) < tol))
