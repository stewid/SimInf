## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2020 Stefan Widgren
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

u0 <- data.frame(S  = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
                 I  = c(0, 0, 0, 0, 0, 0, 0, 0, 0))

## Place nodes in a grid
distance <- expand.grid(x = seq_len(3),
                        y = seq_len(3))
distance <- distance_matrix(distance$x, distance$y, 2)

## Check missing columns in u0
res <- assertError(SISe_sp(u0 = u0[, "I", drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- assertError(SISe_sp(u0 = u0[, "S", drop = FALSE]))
check_error(res, "Missing columns in u0.")

## Check default phi
res <- SISe_sp(u0       = u0,
               tspan    = seq_len(10) - 1,
               events   = NULL,
               upsilon  = 0.0357,
               gamma    = 0.1,
               alpha    = 1.0,
               beta_t1  = 0.19,
               beta_t2  = 0.085,
               beta_t3  = 0.075,
               beta_t4  = 0.185,
               end_t1   = 91,
               end_t2   = 182,
               end_t3   = 273,
               end_t4   = 365,
               coupling = 0.0005,
               distance = distance)
stopifnot(identical(res@v0,
                    structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0),
                              .Dim = c(1L, 9L),
                              .Dimnames = list("phi", NULL))))

## Check missing upsilon
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'upsilon' must be numeric of length 1.")

## Check missing gamma
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'gamma' must be numeric of length 1.")

## Check missing alpha
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'alpha' must be numeric of length 1.")

## Check missing beta_t1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t1' must be numeric of length 1.")

## Check missing beta_t2
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t2' must be numeric of length 1.")

## Check missing beta_t3
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t3' must be numeric of length 1.")

## Check missing beta_t4
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t4' must be numeric of length 1.")

## Check missing end_t1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t1' is missing.")

## Check missing end_t2
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t2' is missing.")

## Check missing end_t3
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t3' is missing.")

## Check missing end_t4
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t4' is missing.")

## Check missing coupling
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           distance = distance))
check_error(res, "'coupling' must be numeric of length 1.")

## Check missing distance
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005))
check_error(res, "'distance' is missing.")

## Check negative distance
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = -distance))
check_error(res, "All values in the 'distance' matrix must be >= 0.")

## Check non-numeric upsilon
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = "0.0357",
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'upsilon' must be numeric of length 1.")

## Check non-numeric gamma
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = "0.1",
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'gamma' must be numeric of length 1.")

## Check non-numeric alpha
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = "1.0",
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'alpha' must be numeric of length 1.")

## Check non-numeric beta_t1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = "0.19",
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t1' must be numeric of length 1.")

## Check non-numeric beta_t2
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = "0.085",
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t2' must be numeric of length 1.")

## Check non-numeric beta_t3
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = "0.075",
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t3' must be numeric of length 1.")

## Check non-numeric beta_t4
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = "0.185",
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t4' must be numeric of length 1.")

## Check non-integer end_t1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = "91",
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t1' must be integer.")

res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91.5,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t1' must be integer.")

## Check non-integer end_t2
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = "182",
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t2' must be integer.")

res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182.5,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t2' must be integer.")

## Check non-integer end_t3
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = "273",
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t3' must be integer.")

res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273.5,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t3' must be integer.")

## Check non-integer end_t4
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = "365",
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t4' must be integer.")

res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365.5,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t4' must be integer.")

## Check non-numeric coupling
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = "0.0005",
                           distance = distance))
check_error(res, "'coupling' must be numeric of length 1.")

## Check that length of upsilon equals 1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = c(0.0357, 0.0357),
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'upsilon' must be numeric of length 1.")

## Check that length of gamma equals 1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = c(0.1, 0.1),
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'gamma' must be numeric of length 1.")

## Check that length of alpha equals 1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = c(1.0, 1.0),
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'alpha' must be numeric of length 1.")

## Check that length of beta_t1 equals 1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = c(0.19, 0.19),
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t1' must be numeric of length 1.")

## Check that length of beta_t2 equals 1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = c(0.085, 0.085),
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t2' must be numeric of length 1.")

## Check that length of beta_t3 equals 1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = c(0.075, 0.075),
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t3' must be numeric of length 1.")

## Check that length of beta_t4 equals 1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = c(0.185, 0.185),
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'beta_t4' must be numeric of length 1.")

## Check that length of coupling equals 1
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = c(0.0005, 0.0005),
                           distance = distance))
check_error(res, "'coupling' must be numeric of length 1.")

## Check interval endpoints
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = -1,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t1' must be greater than or equal to '0'.")

res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 18,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t1' must be less than 'end_t2'.")

res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 173,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(
    res,
    "'end_t2' must be less than 'end_t3' or 'end_t3' less than 'end_t1'.")

res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 365,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t3' must be less than '364'.")

res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = -1,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t4' must be greater than or equal to '0'.")

res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 366,
                           coupling = 0.0005,
                           distance = distance))
check_error(res, "'end_t4' must be less than or equal to '365'.")

res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 4:12,
                           end_t2   = 5:13,
                           end_t3   = c(8:15, 19),
                           end_t4   = c(2, 11:18),
                           coupling = 0.0005,
                           distance = distance))
check_error(
    res,
    "'end_t4' must be less than 'end_t1' or greater than 'end_t3'.")

## Check distance matrix
res <- assertError(SISe_sp(u0       = u0,
                           tspan    = seq_len(10) - 1,
                           events   = NULL,
                           phi      = rep(1, nrow(u0)),
                           upsilon  = 0.0357,
                           gamma    = 0.1,
                           alpha    = 1.0,
                           beta_t1  = 0.19,
                           beta_t2  = 0.085,
                           beta_t3  = 0.075,
                           beta_t4  = 0.185,
                           end_t1   = 91,
                           end_t2   = 182,
                           end_t3   = 273,
                           end_t4   = 365,
                           coupling = 0.0005,
                           distance = as.matrix(distance)))
check_error(res, "The 'distance' argument must be of type 'dgCMatrix'.")

## Check that non-data.frame u0 works
SISe_sp(u0       = cbind(S = 0:8, I = rep(0, 9)),
        tspan    = seq_len(10) - 1,
        events   = NULL,
        phi      = rep(0, 9),
        upsilon  = 0.0357,
        gamma    = 0.1,
        alpha    = 1.0,
        beta_t1  = 0.19,
        beta_t2  = 0.085,
        beta_t3  = 0.075,
        beta_t4  = 0.185,
        end_t1   = 91,
        end_t2   = 182,
        end_t3   = 273,
        end_t4   = 365,
        coupling = 0.0005,
        distance = distance)

## Check extraction of data from 'suscpetible', and 'infected'
## compartments
model <- SISe_sp(u0       = u0,
                 tspan    = seq_len(10) - 1,
                 events   = NULL,
                 phi      = rep(0, nrow(u0)),
                 upsilon  = 0.0357,
                 gamma    = 0.1,
                 alpha    = 1.0,
                 beta_t1  = 0.19,
                 beta_t2  = 0.085,
                 beta_t3  = 0.075,
                 beta_t4  = 0.185,
                 end_t1   = 91,
                 end_t2   = 182,
                 end_t3   = 273,
                 end_t4   = 365,
                 coupling = 0.0005,
                 distance = distance)

result <- run(model)

S_expected <- structure(c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                          0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                          0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                          0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                          0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                          0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                          0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                          0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                          0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                          0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L),
                        .Dim = 9:10)
S_observed <- trajectory(result, compartments = "S", as.is = TRUE)
stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = 9:10)
I_observed <- trajectory(result, compartments = "I", as.is = TRUE)
stopifnot(identical(I_observed, I_expected))

## Check SISe_sp plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SISe_sp run function fails for misspecified SISe_sp model
res <- assertError(.Call(SimInf:::SISe_sp_run, NULL, NULL))
check_error(res, "Invalid model.")

res <- assertError(.Call(SimInf:::SISe_sp_run, "SISe_sp", NULL))
check_error(res, "Invalid model.")

## Check error non-finite v
model <- SISe_sp(u0       = u0,
                 tspan    = seq_len(10) - 1,
                 events   = NULL,
                 phi      = rep(1, nrow(u0)),
                 upsilon  = 0.0357,
                 gamma    = 0.1,
                 alpha    = 1.0,
                 beta_t1  = 0.19,
                 beta_t2  = 0.085,
                 beta_t3  = 0.075,
                 beta_t4  = 0.185,
                 end_t1   = 91,
                 end_t2   = 182,
                 end_t3   = 273,
                 end_t4   = 365,
                 coupling = 0.0005,
                 distance = distance)
model@gdata["beta_t1"] <- Inf
model@gdata["beta_t2"] <- Inf
model@gdata["beta_t3"] <- Inf
model@gdata["beta_t4"] <- Inf
res <- assertError(run(model))
check_error(res, "The continuous state 'v' is not finite.")

## Check negative v
u0 <- data.frame(S  = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
                 I  = c(10, 10, 10, 10, 10, 10, 10, 10, 10))
model <- SISe_sp(u0       = u0,
                 tspan    = seq_len(10) - 1,
                 events   = NULL,
                 phi      = rep(1, nrow(u0)),
                 upsilon  = 0,
                 gamma    = 0,
                 alpha    = -1.0,
                 beta_t1  = 0.19,
                 beta_t2  = 0.085,
                 beta_t3  = 0.075,
                 beta_t4  = 0.185,
                 end_t1   = 91,
                 end_t2   = 182,
                 end_t3   = 273,
                 end_t4   = 365,
                 coupling = 0.0005,
                 distance = distance)
res <- assertError(run(model))
check_error(res, "The continuous state 'v' is negative.")
