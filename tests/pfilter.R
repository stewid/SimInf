## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2021 Stefan Widgren
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

## Specify the number of threads to use.
set_num_threads(1)

## For debugging
sessionInfo()

## Create an SIR model object.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5,
             beta = 0.16,
             gamma = 0.077)

## Check that a missing 'time' column in data raises an error.
res <- assertError(SimInf:::pfilter_tspan(model, data.frame()))
check_error(res, "Missing 'time' column in data.")

## Check that a non-integer value in the 'data$time' column raises an
## error.
res <- assertError(SimInf:::pfilter_tspan(model, data.frame(time = 1.1)))
check_error(res, "'data$time' must be integer.")

## Check that a NA value in the 'data$time' column raises an error.
res <- assertError(SimInf:::pfilter_tspan(model, data.frame(time = NA)))
check_error(res, "'data$time' must be integer.")

## Check that data$time[1] < model@tspan[1] raises an error.
res <- assertError(SimInf:::pfilter_tspan(model, data.frame(time = 0)))
check_error(res, "data$time[1] must be >= tspan[1].")

stopifnot(identical(
    SimInf:::pfilter_tspan(model, data.frame(time = 1:3)),
    list(1L, 1:2, 2:3)))

stopifnot(identical(
    SimInf:::pfilter_tspan(model, data.frame(time = 2:3)),
    list(c(1, 2), 2:3)))
