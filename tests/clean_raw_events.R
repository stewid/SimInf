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

## Check to pass vectors of different lengths.
res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    integer(0),
    c(0L, 3L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'id' must be an integer vector with length > 0.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    c(1L, 1L),
    c(0L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'event' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    c(1L, 1L),
    c(0L, 3L),
    c(1L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'time' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    c(1L, 1L),
    c(0L, 3L),
    c(1L, 2L),
    c(1L),
    c(0L, 2L)))
check_error(res, "'node' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    c(1L, 1L),
    c(0L, 3L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L)))
check_error(res, "'dest' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    c(1L, 1L),
    c(0L, 2L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'event[2]' is invalid.")
