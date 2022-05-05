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

## Specify the number of threads to use.
set_num_threads(1)

## For debugging
sessionInfo()

## Check invalid u0
res <- assertError(SIS(u0 = "u0"))
check_error(res, "Missing columns in u0.")

u0 <- data.frame(S  = c(0, 1, 2, 3, 4, 5),
                 I  = c(0, 0, 0, 0, 0, 0))

## Check missing columns in u0
res <- assertError(SIS(u0 = u0[, "I", drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- assertError(SIS(u0 = u0[, "S", drop = FALSE]))
check_error(res, "Missing columns in u0.")

## Check missing beta
res <- assertError(SIS(u0     = u0,
                       tspan  = seq_len(10) - 1,
                       events = NULL,
                       gamma  = 0.0357))
check_error(res, "'beta' must be numeric of length 1 or 'nrow(u0)'.")

## Check non-numeric beta
res <- assertError(SIS(u0     = u0,
                       tspan  = seq_len(10) - 1,
                       events = NULL,
                       beta   = "0.1",
                       gamma  = 0.0357))
check_error(res, "'beta' must be numeric of length 1 or 'nrow(u0)'.")

## Check length of beta
res <- assertError(SIS(u0     = u0,
                       tspan  = seq_len(10) - 1,
                       events = NULL,
                       beta   = rep(0.1, nrow(u0) + 1),
                       gamma  = 0.0357))
check_error(res, "'beta' must be numeric of length 1 or 'nrow(u0)'.")

## Check missing gamma
res <- assertError(SIS(u0     = u0,
                       tspan  = seq_len(10) - 1,
                       events = NULL,
                       beta   = 0.0357))
check_error(res, "'gamma' must be numeric of length 1 or 'nrow(u0)'.")

## Check non-numeric gamma
res <- assertError(SIS(u0     = u0,
                       tspan  = seq_len(10) - 1,
                       events = NULL,
                       beta   = 0.0357,
                       gamma  = "0.1"))
check_error(res, "'gamma' must be numeric of length 1 or 'nrow(u0)'.")

## Check length of gamma
res <- assertError(SIS(u0     = u0,
                       tspan  = seq_len(10) - 1,
                       events = NULL,
                       beta   = 0.0357,
                       gamma  = rep(0.1, nrow(u0) + 1)))
check_error(res, "'gamma' must be numeric of length 1 or 'nrow(u0)'.")
