## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2019 Stefan Widgren
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

library("SimInf")
source("util/check.R")

## Specify the number of threads to use.
set_num_threads(1)

## For debugging
sessionInfo()

## Check invalid u0
res <- tools::assertError(SISe3_sp(u0 = "u0"))
check_error(res, "Missing columns in u0.")

u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
                 I_1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
                 I_2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
                 I_3 = c(0, 0, 0, 0, 0, 0, 0, 0, 0))

## Place nodes in a grid
distance <- expand.grid(x = seq_len(3),
                        y = seq_len(3))
distance <- distance_matrix(distance$x, distance$y, 2)

## Check missing columns in u0
res <- tools::assertError(
    SISe3_sp(u0 = u0[, c("I_1", "S_2", "I_2", "S_3", "I_3")]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(
    SISe3_sp(u0 = u0[, c("S_1", "S_2", "I_2", "S_3", "I_3")]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(
    SISe3_sp(u0 = u0[, c("S_1", "I_1", "I_2", "S_3", "I_3")]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(
    SISe3_sp(u0 = u0[, c("S_1", "I_1", "S_2", "S_3", "I_3")]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(
    SISe3_sp(u0 = u0[, c("S_1", "I_1", "S_2", "I_2", "I_3")]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(
    SISe3_sp(u0 = u0[, c("S_1", "I_1", "S_2", "I_2", "S_3")]))
check_error(res, "Missing columns in u0.")

## Check default phi
res <- SISe3_sp(u0        = u0,
                tspan     = seq_len(10) - 1,
                events    = NULL,
                upsilon_1 = 0.0357,
                upsilon_2 = 0.0357,
                upsilon_3 = 0.00935,
                gamma_1   = 0.1,
                gamma_2   = 0.1,
                gamma_3   = 0.1,
                alpha     = 1.0,
                beta_t1   = 0.19,
                beta_t2   = 0.085,
                beta_t3   = 0.075,
                beta_t4   = 0.185,
                end_t1    = 91,
                end_t2    = 182,
                end_t3    = 273,
                end_t4    = 365,
                coupling  = 0.0005,
                distance  = distance)
stopifnot(identical(res@v0,
                    structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0),
                              .Dim = c(1L, 9L),
                              .Dimnames = list("phi", NULL))))

## Check missing upsilon_1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'upsilon_1' is missing.")

## Check missing upsilon_2
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'upsilon_2' is missing.")

## Check missing upsilon_3
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'upsilon_3' is missing.")

## Check missing gamma_1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'gamma_1' is missing.")

## Check missing gamma_2
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'gamma_2' is missing.")

## Check missing gamma_3
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'gamma_3' is missing.")

## Check missing alpha
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'alpha' is missing.")

## Check missing beta_t1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t1' is missing.")

## Check missing beta_t2
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t2' is missing.")

## Check missing beta_t3
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t3' is missing.")

## Check missing beta_t4
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t4' is missing.")

## Check missing end_t1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t1' is missing.")

## Check missing end_t2
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t2' is missing.")

## Check missing end_t3
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t3' is missing.")

## Check missing end_t4
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t4' is missing.")

## Check missing coupling
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   distance  = distance))
check_error(res, "'coupling' is missing.")

## Check missing distance
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005))
check_error(res, "'distance' is missing.")

## Check negative distance
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = -distance))
check_error(res, "All values in the 'distance' matrix must be >= 0.")

## Check non-numeric upsilon_1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = "0.0357",
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'upsilon_1' must be numeric.")

## Check non-numeric upsilon_2
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = "0.0357",
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'upsilon_2' must be numeric.")

## Check non-numeric upsilon_3
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = "0.00935",
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'upsilon_3' must be numeric.")

## Check non-numeric gamma_1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = "0.1",
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'gamma_1' must be numeric.")

## Check non-numeric gamma_2
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = "0.1",
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'gamma_2' must be numeric.")

## Check non-numeric gamma_3
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = "0.1",
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'gamma_3' must be numeric.")

## Check non-numeric alpha
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = "1.0",
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'alpha' must be numeric.")

## Check non-numeric beta_t1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = "0.19",
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t1' must be numeric.")

## Check non-numeric beta_t2
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = "0.085",
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t2' must be numeric.")

## Check non-numeric beta_t3
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = "0.075",
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t3' must be numeric.")

## Check non-numeric beta_t4
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = "0.185",
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t4' must be numeric.")

## Check non-integer end_t1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = "91",
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t1' must be integer.")

res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91.5,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t1' must be integer.")

## Check non-integer end_t2
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = "182",
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t2' must be integer.")

res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182.5,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t2' must be integer.")

## Check non-integer end_t3
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = "273",
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t3' must be integer.")

res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273.5,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t3' must be integer.")

## Check non-integer end_t4
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = "365",
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t4' must be integer.")

res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365.5,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t4' must be integer.")

## Check non-numeric coupling
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = "0.0005",
                                   distance  = distance))
check_error(res, "'coupling' must be numeric.")

## Check that length of upsilon_1 equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = c(0.0357, 0.0357),
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'upsilon_1' must be of length 1.")

## Check that length of upsilon_2 equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = c(0.0357, 0.0357),
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'upsilon_2' must be of length 1.")

## Check that length of upsilon_3 equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = c(0.00935, 0.00935),
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'upsilon_3' must be of length 1.")

## Check that length of gamma_1 equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = c(0.1, 0.1),
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'gamma_1' must be of length 1.")

## Check that length of gamma_2 equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = c(0.1, 0.1),
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'gamma_2' must be of length 1.")

## Check that length of gamma_3 equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = c(0.1, 0.1),
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'gamma_3' must be of length 1.")

## Check that length of alpha equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = c(1.0, 1.0),
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'alpha' must be of length 1.")

## Check that length of beta_t1 equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = c(0.19, 0.19),
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t1' must be of length 1.")

## Check that length of beta_t2 equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = c(0.085, 0.085),
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t2' must be of length 1.")

## Check that length of beta_t3 equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = c(0.075, 0.075),
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t3' must be of length 1.")

## Check that length of beta_t4 equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = c(0.185, 0.185),
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'beta_t4' must be of length 1.")

## Check that length of coupling equals 1
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = c(0.0005, 0.0005),
                                   distance  = distance))
check_error(res, "'coupling' must be of length 1.")

## Check interval endpoints
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = -1,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t1' must be greater than or equal to '0'.")

res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 18,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t1' must be less than 'end_t2'.")

res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 173,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(
    res,
    "'end_t2' must be less than 'end_t3' or 'end_t3' less than 'end_t1'.")

res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 365,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t3' must be less than '364'.")

res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = -1,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t4' must be greater than or equal to '0'.")

res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 366,
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(res, "'end_t4' must be less than or equal to '365'.")

res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 4:12,
                                   end_t2    = 5:13,
                                   end_t3    = c(8:15, 19),
                                   end_t4    = c(2, 11:18),
                                   coupling  = 0.0005,
                                   distance  = distance))
check_error(
    res,
    "'end_t4' must be less than 'end_t1' or greater than 'end_t3'.")

## Check distance matrix
res <- tools::assertError(SISe3_sp(u0        = u0,
                                   tspan     = seq_len(10) - 1,
                                   events    = NULL,
                                   phi       = rep(1, nrow(u0)),
                                   upsilon_1 = 0.0357,
                                   upsilon_2 = 0.0357,
                                   upsilon_3 = 0.00935,
                                   gamma_1   = 0.1,
                                   gamma_2   = 0.1,
                                   gamma_3   = 0.1,
                                   alpha     = 1.0,
                                   beta_t1   = 0.19,
                                   beta_t2   = 0.085,
                                   beta_t3   = 0.075,
                                   beta_t4   = 0.185,
                                   end_t1    = 91,
                                   end_t2    = 182,
                                   end_t3    = 273,
                                   end_t4    = 365,
                                   coupling  = 0.0005,
                                   distance  = as.matrix(distance)))
check_error(res, "The 'distance' argument must be of type 'dgCMatrix'.")

## Check extraction of data from 'suscpetible', and 'infected'
## compartments
model <- SISe3_sp(u0        = u0,
                  tspan     = seq_len(10) - 1,
                  events    = NULL,
                  phi       = rep(0, nrow(u0)),
                  upsilon_1 = 0.0357,
                  upsilon_2 = 0.0357,
                  upsilon_3 = 0.00935,
                  gamma_1   = 0.1,
                  gamma_2   = 0.1,
                  gamma_3   = 0.1,
                  alpha     = 1.0,
                  beta_t1   = 0.19,
                  beta_t2   = 0.085,
                  beta_t3   = 0.075,
                  beta_t4   = 0.185,
                  end_t1    = 91,
                  end_t2    = 182,
                  end_t3    = 273,
                  end_t4    = 365,
                  coupling  = 0.0005,
                  distance  = distance)

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
S_observed <- trajectory(result, compartments = "S_1", as.is = TRUE)
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
I_observed <- trajectory(result, compartments = "I_1", as.is = TRUE)
stopifnot(identical(I_observed, I_expected))

## Check SISe3_sp plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SISe3_sp run function fails for misspecified SISe3_sp model
res <- tools::assertError(.Call(SimInf:::SISe3_sp_run, NULL, NULL, NULL))
check_error(res, "Invalid model.")

res <- tools::assertError(.Call(SimInf:::SISe3_sp_run, "SISe3_sp", NULL, NULL))
check_error(res, "Invalid model.")

## Check error non-finite v
model <- SISe3_sp(u0        = u0,
                  tspan     = seq_len(10) - 1,
                  events    = NULL,
                  phi       = rep(1, nrow(u0)),
                  upsilon_1 = 0.0357,
                  upsilon_2 = 0.0357,
                  upsilon_3 = 0.00935,
                  gamma_1   = 0.1,
                  gamma_2   = 0.1,
                  gamma_3   = 0.1,
                  alpha     = 1.0,
                  beta_t1   = 0.19,
                  beta_t2   = 0.085,
                  beta_t3   = 0.075,
                  beta_t4   = 0.185,
                  end_t1    = 91,
                  end_t2    = 182,
                  end_t3    = 273,
                  end_t4    = 365,
                  coupling  = 0.0005,
                  distance  = distance)
model@gdata["beta_t1"] <- Inf
model@gdata["beta_t2"] <- Inf
model@gdata["beta_t3"] <- Inf
model@gdata["beta_t4"] <- Inf
res <- tools::assertError(run(model))
check_error(res, "The continuous state 'v' is not finite.")

## Check negative v
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
                 I_1 = c(10, 10, 10, 10, 10, 10, 10, 10, 10),
                 S_2 = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
                 I_2 = c(10, 10, 10, 10, 10, 10, 10, 10, 10),
                 S_3 = c(0, 1, 2, 3, 4, 5, 6, 7, 8),
                 I_3 = c(10, 10, 10, 10, 10, 10, 10, 10, 10))
model <- SISe3_sp(u0        = u0,
                  tspan     = seq_len(10) - 1,
                  events    = NULL,
                  phi       = rep(1, nrow(u0)),
                  upsilon_1 = 0,
                  upsilon_2 = 0,
                  upsilon_3 = 0,
                  gamma_1   = 0,
                  gamma_2   = 0,
                  gamma_3   = 0,
                  alpha     = -1.0,
                  beta_t1   = 0.19,
                  beta_t2   = 0.085,
                  beta_t3   = 0.075,
                  beta_t4   = 0.185,
                  end_t1    = 91,
                  end_t2    = 182,
                  end_t3    = 273,
                  end_t4    = 365,
                  coupling  = 0.0005,
                  distance  = distance)
res <- tools::assertError(run(model))
check_error(res, "The continuous state 'v' is negative.")
