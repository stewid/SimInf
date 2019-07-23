## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015 - 2019  Stefan Engblom
## Copyright (C) 2015 - 2019  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
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
res <- tools::assertError(SISe3(u0 = "u0"))
check_error(res, "Missing columns in u0.")

u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

## Check missing columns in u0
res <- tools::assertError(
    SISe3(u0 = u0[, c("I_1", "S_2", "I_2", "S_3", "I_3")]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(
    SISe3(u0 = u0[, c("S_1", "S_2", "I_2", "S_3", "I_3")]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(
    SISe3(u0 = u0[, c("S_1", "I_1", "I_2", "S_3", "I_3")]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(
    SISe3(u0 = u0[, c("S_1", "I_1", "S_2", "S_3", "I_3")]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(
    SISe3(u0 = u0[, c("S_1", "I_1", "S_2", "I_2", "I_3")]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(
    SISe3(u0 = u0[, c("S_1", "I_1", "S_2", "I_2", "S_3")]))
check_error(res, "Missing columns in u0.")

## Check default phi
res <- SISe3(u0        = u0,
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
             epsilon   = 0.000011)
stopifnot(identical(res@v0,
                    structure(c(0, 0, 0, 0, 0, 0),
                              .Dim = c(1L, 6L),
                              .Dimnames = list("phi", NULL))))

## Check missing upsilon_1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'upsilon_1' is missing.")

## Check missing upsilon_2
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'upsilon_2' is missing.")

## Check missing upsilon_3
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'upsilon_3' is missing.")

## Check missing gamma_1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'gamma_1' is missing.")

## Check missing gamma_2
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'gamma_2' is missing.")

## Check missing gamma_3
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'gamma_3' is missing.")

## Check missing alpha
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'alpha' is missing.")

## Check missing beta_t1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t1' is missing.")

## Check missing beta_t2
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t2' is missing.")

## Check missing beta_t3
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t3' is missing.")

## Check missing beta_t4
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t4' is missing.")

## Check missing end_t1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t1' is missing.")

## Check missing end_t2
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t2' is missing.")

## Check missing end_t3
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t3' is missing.")

## Check missing end_t4
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t4' is missing.")

## Check missing epsilon
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                end_t4    = 365))
check_error(res, "'epsilon' is missing.")

## Check non-numeric upsilon_1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'upsilon_1' must be numeric.")

## Check non-numeric upsilon_2
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'upsilon_2' must be numeric.")

## Check non-numeric upsilon_3
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'upsilon_3' must be numeric.")

## Check non-numeric gamma_1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'gamma_1' must be numeric.")

## Check non-numeric gamma_2
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'gamma_2' must be numeric.")

## Check non-numeric gamma_3
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'gamma_3' must be numeric.")

## Check non-numeric alpha
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'alpha' must be numeric.")

## Check non-numeric beta_t1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t1' must be numeric.")

## Check non-numeric beta_t2
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t2' must be numeric.")

## Check non-numeric beta_t3
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t3' must be numeric.")

## Check non-numeric beta_t4
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t4' must be numeric.")

## Check non-integer end_t1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t1' must be integer.")

res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t1' must be integer.")

## Check non-integer end_t2
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t2' must be integer.")

res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t2' must be integer.")

## Check non-integer end_t3
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t3' must be integer.")

res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t3' must be integer.")

## Check non-integer end_t4
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t4' must be integer.")

res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t4' must be integer.")

## Check non-numeric epsilon
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = "0.000011"))
check_error(res, "'epsilon' must be numeric.")

## Check that length of upsilon_1 equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'upsilon_1' must be of length 1.")

## Check that length of upsilon_2 equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'upsilon_2' must be of length 1.")

## Check that length of upsilon_3 equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'upsilon_3' must be of length 1.")

## Check that length of gamma_1 equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'gamma_1' must be of length 1.")

## Check that length of gamma_2 equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'gamma_2' must be of length 1.")

## Check that length of gamma_3 equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'gamma_3' must be of length 1.")

## Check that length of alpha equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'alpha' must be of length 1.")

## Check that length of beta_t1 equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t1' must be of length 1.")

## Check that length of beta_t2 equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t2' must be of length 1.")

## Check that length of beta_t3 equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t3' must be of length 1.")

## Check that length of beta_t4 equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'beta_t4' must be of length 1.")

## Check that length of epsilon equals 1
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = c(0.000011, 0.000011)))
check_error(res, "'epsilon' must be of length 1.")

## Check interval endpoints
res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t1' must be greater than or equal to '0'.")

res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t1' must be less than 'end_t2'.")

res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t2' must be less than 'end_t3' or 'end_t3' less than 'end_t1'.")

res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t3' must be less than '364'.")

res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t4' must be greater than or equal to '0'.")

res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                epsilon   = 0.000011))
check_error(res, "'end_t4' must be less than or equal to '365'.")

res <- tools::assertError(SISe3(u0        = u0,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
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
                                end_t1    = 4:9,
                                end_t2    = 5:10,
                                end_t3    = c(8:12, 16),
                                end_t4    = c(2, 11:15),
                                epsilon   = 0.000011))
check_error(res, "'end_t4' must be less than 'end_t1' or greater than 'end_t3'.")

## Check extraction of data from 'suscpetible', and 'infected'
## compartments
model <- SISe3(u0        = u0,
               tspan     = seq_len(10) - 1,
               events    = NULL,
               phi       = rep(0, 6),
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
               epsilon   = 0.000011)

set.seed(123)
result <- run(model)

S_expected <- structure(c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L,
                          1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
                          2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L,
                          3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
                          4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                        .Dim = c(6L, 10L))
S_observed <- trajectory(result, compartments = "S_1", as.is = TRUE)
stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L))
I_observed <- trajectory(result, compartments = "I_1", as.is = TRUE)
stopifnot(identical(I_observed, I_expected))

## Check SISe3 plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SISe3 run function fails for misspecified SISe3 model
res <- tools::assertError(.Call("SISe3_run", NULL, NULL, NULL, PACKAGE = "SimInf"))
check_error(res, "Invalid model.")

res <- tools::assertError(.Call("SISe3_run", "SISe3", NULL, NULL, PACKAGE = "SimInf"))
check_error(res, "Invalid model.")

## Check error non-finite v
model <- SISe3(u0        = u0,
               tspan     = seq_len(10) - 1,
               events    = NULL,
               phi       = rep(1, 6),
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
               epsilon   = 0.000011)
model@gdata["beta_t1"] <- Inf
model@gdata["beta_t2"] <- Inf
model@gdata["beta_t3"] <- Inf
model@gdata["beta_t4"] <- Inf
res <- tools::assertError(run(model))
check_error(res, "The continuous state 'v' is not finite.")

## Check negative v
model <- SISe3(u0        = u0,
               tspan     = seq_len(10) - 1,
               events    = NULL,
               phi       = rep(1, 6),
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
               epsilon   = -10.000011)
res <- tools::assertError(run(model))
check_error(res, "The continuous state 'v' is negative.")
