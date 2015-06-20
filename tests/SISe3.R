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
tools::assertError(SISe3(init = init[, c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3")]))
tools::assertError(SISe3(init = init[, c("id", "I_1", "S_2", "I_2", "S_3", "I_3")]))
tools::assertError(SISe3(init = init[, c("id", "S_1", "S_2", "I_2", "S_3", "I_3")]))
tools::assertError(SISe3(init = init[, c("id", "S_1", "I_I", "I_2", "S_3", "I_3")]))
tools::assertError(SISe3(init = init[, c("id", "S_1", "I_I", "S_2", "S_3", "I_3")]))
tools::assertError(SISe3(init = init[, c("id", "S_1", "I_I", "S_2", "I_2", "I_3")]))
tools::assertError(SISe3(init = init[, c("id", "S_1", "I_I", "S_2", "I_2", "S_3")]))

## Check missing upsilon_1
tools::assertError(SISe3(init      = init,
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

## Check missing upsilon_2
tools::assertError(SISe3(init      = init,
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

## Check missing upsilon_3
tools::assertError(SISe3(init      = init,
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

## Check missing gamma_1
tools::assertError(SISe3(init      = init,
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

## Check missing gamma_2
tools::assertError(SISe3(init      = init,
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

## Check missing gamma_3
tools::assertError(SISe3(init      = init,
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

## Check missing alpha
tools::assertError(SISe3(init      = init,
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

## Check missing beta_q1
tools::assertError(SISe3(init      = init,
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

## Check missing beta_q2
tools::assertError(SISe3(init      = init,
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

## Check missing beta_q3
tools::assertError(SISe3(init      = init,
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

## Check missing beta_q4
tools::assertError(SISe3(init      = init,
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

## Check missing epsilon
tools::assertError(SISe3(init      = init,
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

## Check non-numeric upsilon_1
tools::assertError(SISe3(init      = init,
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

## Check non-numeric upsilon_2
tools::assertError(SISe3(init      = init,
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

## Check non-numeric upsilon_3
tools::assertError(SISe3(init      = init,
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

## Check non-numeric gamma_1
tools::assertError(SISe3(init      = init,
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

## Check non-numeric gamma_2
tools::assertError(SISe3(init      = init,
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

## Check non-numeric gamma_3
tools::assertError(SISe3(init      = init,
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

## Check non-numeric alpha
tools::assertError(SISe3(init      = init,
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

## Check non-numeric beta_q1
tools::assertError(SISe3(init      = init,
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

## Check non-numeric beta_q2
tools::assertError(SISe3(init      = init,
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

## Check non-numeric beta_q3
tools::assertError(SISe3(init      = init,
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

## Check non-numeric beta_q4
tools::assertError(SISe3(init      = init,
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

## Check non-numeric epsilon
tools::assertError(SISe3(init      = init,
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

## Check that length of upsilon_1 equals 1
tools::assertError(SISe3(init      = init,
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

## Check that length of upsilon_2 equals 1
tools::assertError(SISe3(init      = init,
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

## Check that length of upsilon_3 equals 1
tools::assertError(SISe3(init      = init,
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

## Check that length of gamma_1 equals 1
tools::assertError(SISe3(init      = init,
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

## Check that length of gamma_2 equals 1
tools::assertError(SISe3(init      = init,
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

## Check that length of gamma_3 equals 1
tools::assertError(SISe3(init      = init,
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

## Check that length of alpha equals 1
tools::assertError(SISe3(init      = init,
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

## Check that length of epsilon equals 1
tools::assertError(SISe3(init      = init,
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

## Check that length of beta_q1 equals 1 or nrow(init)
tools::assertError(SISe3(init      = init,
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
                         beta_q1   = rep(0.19, nrow(init) + 1),
                         beta_q2   = 0.085,
                         beta_q3   = 0.075,
                         beta_q4   = 0.185,
                         epsilon   = 0.000011))

## Check that length of beta_q2 equals 1 or nrow(init)
tools::assertError(SISe3(init      = init,
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
                         beta_q2   = rep(0.085, nrow(init) + 1),
                         beta_q3   = 0.075,
                         beta_q4   = 0.185,
                         epsilon   = 0.000011))

## Check that length of beta_q3 equals 1 or nrow(init)
tools::assertError(SISe3(init      = init,
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
                         beta_q3   = rep(0.075, nrow(init) + 1),
                         beta_q4   = 0.185,
                         epsilon   = 0.000011))

## Check that length of beta_q4 equals 1 or nrow(init)
tools::assertError(SISe3(init      = init,
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
                         beta_q4   = rep(0.185, nrow(init) + 1),
                         epsilon   = 0.000011))

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
