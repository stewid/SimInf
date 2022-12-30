## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2022 Stefan Widgren
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
library(tools)
source("util/check.R")

res <- assertError(SimInf:::KLIEP(xnu = 1:5))
check_error(res, "'xnu' must be a numeric matrix.")

res <- assertError(SimInf:::KLIEP(xnu = matrix(letters)))
check_error(res, "'xnu' must be a numeric matrix.")

res <- assertError(SimInf:::KLIEP(xnu = matrix(1:5), xde = 1:5))
check_error(res, "'xde' must be a numeric matrix.")

res <- assertError(SimInf:::KLIEP(xnu = matrix(1:5), xde = matrix(letters)))
check_error(res, "'xde' must be a numeric matrix.")

res <- assertError(SimInf:::KLIEP(xnu = matrix(1:5), xde = matrix(1:10, ncol = 2)))
check_error(res, "'xnu' and 'xde' must have the same dimension.")
