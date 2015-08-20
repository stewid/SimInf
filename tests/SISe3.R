## siminf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015  Stefan Engblom
## Copyright (C) 2015  Stefan Widgren
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
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(siminf)

init <- structure(list(id  = c(0, 1, 2, 3, 4, 5),
                       S_1 = c(0, 1, 2, 3, 4, 5),
                       I_1 = c(0, 0, 0, 0, 0, 0),
                       S_2 = c(0, 1, 2, 3, 4, 5),
                       I_2 = c(0, 0, 0, 0, 0, 0),
                       S_3 = c(0, 1, 2, 3, 4, 5),
                       I_3 = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -6L), class = "data.frame")

## Check missing columns in init
res <- tools::assertError(
    SISe3(init = init[, c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3")]))
stopifnot(length(grep("Missing columns in init",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    SISe3(init = init[, c("id", "I_1", "S_2", "I_2", "S_3", "I_3")]))
stopifnot(length(grep("Missing columns in init",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    SISe3(init = init[, c("id", "S_1", "S_2", "I_2", "S_3", "I_3")]))
stopifnot(length(grep("Missing columns in init",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    SISe3(init = init[, c("id", "S_1", "I_1", "I_2", "S_3", "I_3")]))
stopifnot(length(grep("Missing columns in init",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    SISe3(init = init[, c("id", "S_1", "I_1", "S_2", "S_3", "I_3")]))
stopifnot(length(grep("Missing columns in init",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    SISe3(init = init[, c("id", "S_1", "I_1", "S_2", "I_2", "I_3")]))
stopifnot(length(grep("Missing columns in init",
                      res[[1]]$message)) > 0)

res <- tools::assertError(
    SISe3(init = init[, c("id", "S_1", "I_1", "S_2", "I_2", "S_3")]))
stopifnot(length(grep("Missing columns in init",
                      res[[1]]$message)) > 0)

## Check missing upsilon_1
res <- tools::assertError(SISe3(init      = init,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
                                upsilon_2 = 0.0357,
                                upsilon_3 = 0.00935,
                                gamma_1   = 0.1,
                                gamma_2   = 0.1,
                                gamma_3   = 0.1,
                                alpha     = 1.0,
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'upsilon_1' is missing",
                      res[[1]]$message)) > 0)

## Check missing upsilon_2
res <- tools::assertError(SISe3(init      = init,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
                                upsilon_1 = 0.0357,
                                upsilon_3 = 0.00935,
                                gamma_1   = 0.1,
                                gamma_2   = 0.1,
                                gamma_3   = 0.1,
                                alpha     = 1.0,
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'upsilon_2' is missing",
                      res[[1]]$message)) > 0)

## Check missing upsilon_3
res <- tools::assertError(SISe3(init      = init,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
                                upsilon_1 = 0.0357,
                                upsilon_2 = 0.0357,
                                gamma_1   = 0.1,
                                gamma_2   = 0.1,
                                gamma_3   = 0.1,
                                alpha     = 1.0,
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'upsilon_3' is missing",
                      res[[1]]$message)) > 0)

## Check missing gamma_1
res <- tools::assertError(SISe3(init      = init,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
                                upsilon_1 = 0.0357,
                                upsilon_2 = 0.0357,
                                upsilon_3 = 0.00935,
                                gamma_2   = 0.1,
                                gamma_3   = 0.1,
                                alpha     = 1.0,
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'gamma_1' is missing",
                      res[[1]]$message)) > 0)

## Check missing gamma_2
res <- tools::assertError(SISe3(init      = init,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
                                upsilon_1 = 0.0357,
                                upsilon_2 = 0.0357,
                                upsilon_3 = 0.00935,
                                gamma_1   = 0.1,
                                gamma_3   = 0.1,
                                alpha     = 1.0,
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'gamma_2' is missing",
                      res[[1]]$message)) > 0)

## Check missing gamma_3
res <- tools::assertError(SISe3(init      = init,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
                                upsilon_1 = 0.0357,
                                upsilon_2 = 0.0357,
                                upsilon_3 = 0.00935,
                                gamma_1   = 0.1,
                                gamma_2   = 0.1,
                                alpha     = 1.0,
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'gamma_3' is missing",
                      res[[1]]$message)) > 0)

## Check missing alpha
res <- tools::assertError(SISe3(init      = init,
                                tspan     = seq_len(10) - 1,
                                events    = NULL,
                                phi       = rep(1, 6),
                                upsilon_1 = 0.0357,
                                upsilon_2 = 0.0357,
                                upsilon_3 = 0.00935,
                                gamma_1   = 0.1,
                                gamma_2   = 0.1,
                                gamma_3   = 0.1,
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'alpha' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_q1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q1' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_q2
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q2' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_q3
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q3' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_q4
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q4' is missing",
                      res[[1]]$message)) > 0)

## Check missing epsilon
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185))
stopifnot(length(grep("'epsilon' is missing",
                      res[[1]]$message)) > 0)

## Check non-numeric upsilon_1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'upsilon_1' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric upsilon_2
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'upsilon_2' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric upsilon_3
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'upsilon_3' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric gamma_1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'gamma_1' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric gamma_2
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'gamma_2' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric gamma_3
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'gamma_3' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric alpha
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'alpha' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_q1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = "0.19",
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q1' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_q2
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = "0.085",
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q2' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_q3
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = "0.075",
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q3' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_q4
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = "0.185",
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q4' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric epsilon
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = "0.000011"))
stopifnot(length(grep("'epsilon' must be numeric",
                      res[[1]]$message)) > 0)

## Check that length of upsilon_1 equals 1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'upsilon_1' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of upsilon_2 equals 1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'upsilon_2' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of upsilon_3 equals 1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'upsilon_3' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of gamma_1 equals 1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'gamma_1' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of gamma_2 equals 1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'gamma_2' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of gamma_3 equals 1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'gamma_3' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of alpha equals 1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'alpha' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_q1 equals 1 or nrow(init)
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = c(0.19, 0.19),
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q1' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_q2 equals 1 or nrow(init)
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = c(0.085, 0.085),
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q2' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_q3 equals 1 or nrow(init)
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = c(0.075, 0.075),
                                beta_q4   = 0.185,
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q3' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_q4 equals 1 or nrow(init)
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = c(0.185, 0.185),
                                epsilon   = 0.000011))
stopifnot(length(grep("'beta_q4' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of epsilon equals 1
res <- tools::assertError(SISe3(init      = init,
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
                                beta_q1   = 0.19,
                                beta_q2   = 0.085,
                                beta_q3   = 0.075,
                                beta_q4   = 0.185,
                                epsilon   = c(0.000011, 0.000011)))
stopifnot(length(grep("'epsilon' must be of length 1",
                      res[[1]]$message)) > 0)

## Check 'suscpetible' and 'infected' methods
model <- SISe3(init      = init,
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
               beta_q1   = 0.19,
               beta_q2   = 0.085,
               beta_q3   = 0.075,
               beta_q4   = 0.185,
               epsilon   = 0.000011)

result <- run(model)

S_expected <- structure(c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L,
                          1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
                          2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L,
                          3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
                          4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                        .Dim = c(6L, 10L))

S_observed <- susceptible(result)

stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L))

I_observed <- infected(result)

stopifnot(identical(I_observed, I_expected))

## Check SISe3 plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result, t0 = "2015-01-01")
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SISe3 run function fails for misspecified SISe3 model
res <- tools::assertError(.Call(siminf:::SISe3_run, NULL, NULL, NULL))
stopifnot(length(grep("Invalid SISe3 model",
                      res[[1]]$message)) > 0)
res <- tools::assertError(.Call(siminf:::SISe3_run, "SISe3", NULL, NULL))
stopifnot(length(grep("Invalid SISe3 model",
                      res[[1]]$message)) > 0)

setClass("DummySISe3", slots = c(a = "character"))
model <- new("DummySISe3", a = "SISe3")
res <- tools::assertError(.Call(siminf:::SISe3_run, model, NULL, NULL))
stopifnot(length(grep("Invalid SISe3 model: DummySISe3",
                      res[[1]]$message)) > 0)
