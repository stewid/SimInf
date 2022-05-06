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

## Check running a trajectory
trajectory_exp <- data.frame(
    node = rep(1L, 100),
    time = 1:100,
    S = c(99L, 99L, 98L, 98L, 98L, 98L, 98L, 98L, 98L, 97L, 96L, 97L,
          97L, 97L, 94L, 97L, 97L, 97L, 97L, 98L, 98L, 98L, 97L, 97L,
          97L, 97L, 97L, 96L, 96L, 97L, 96L, 96L, 94L, 94L, 94L, 93L,
          93L, 92L, 93L, 92L, 94L, 94L, 93L, 92L, 93L, 95L, 93L, 93L,
          93L, 93L, 92L, 91L, 91L, 92L, 90L, 89L, 91L, 92L, 92L, 92L,
          92L, 92L, 92L, 89L, 85L, 83L, 81L, 80L, 80L, 79L, 82L, 82L,
          82L, 79L, 77L, 75L, 77L, 77L, 77L, 77L, 77L, 76L, 76L, 76L,
          74L, 74L, 73L, 74L, 73L, 71L, 69L, 72L, 70L, 70L, 73L, 74L,
          76L, 72L, 74L, 71L),
    I = c(1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 4L, 3L, 3L, 3L, 6L,
          3L, 3L, 3L, 3L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 3L,
          4L, 4L, 6L, 6L, 6L, 7L, 7L, 8L, 7L, 8L, 6L, 6L, 7L, 8L, 7L,
          5L, 7L, 7L, 7L, 7L, 8L, 9L, 9L, 8L, 10L, 11L, 9L, 8L, 8L,
          8L, 8L, 8L, 8L, 11L, 15L, 17L, 19L, 20L, 20L, 21L, 18L, 18L,
          18L, 21L, 23L, 25L, 23L, 23L, 23L, 23L, 23L, 24L, 24L, 24L,
          26L, 26L, 27L, 26L, 27L, 29L, 31L, 28L, 30L, 30L, 27L, 26L,
          24L, 28L, 26L, 29L))

model <- SIS(u0 = data.frame(S = 99, I = 1),
             tspan = 1:100,
             events = NULL,
             beta = 0.16,
             gamma = 0.077)

set.seed(22)
trajectory_obs <- trajectory(run(model))
stopifnot(identical(trajectory_obs, trajectory_exp))

## Check data
stopifnot(identical(events_SIS(), events_SISe()))
stopifnot(identical(u0_SIS(), u0_SISe()))
