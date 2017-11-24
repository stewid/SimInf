## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015 - 2017  Stefan Engblom
## Copyright (C) 2015 - 2017  Stefan Widgren
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

library("SimInf")

## For debugging
sessionInfo()

## Check invalid u0
res <- tools::assertError(SISe(u0 = "u0"))
stopifnot(length(grep("'u0' must be a data.frame",
                      res[[1]]$message)) > 0)

u0 <- structure(list(S  = c(0, 1, 2, 3, 4, 5),
                     I  = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S", "I"),
                row.names = c(NA, -6L), class = "data.frame")

## Check missing columns in u0
res <- tools::assertError(SISe(u0 = u0[, "I", drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SISe(u0 = u0[, "S", drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)

## Check phi
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep("a", nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("Invalid 'phi': must be numeric vector",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = matrix(rep(1, nrow(u0))),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("Invalid 'phi': must be numeric vector",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0) - 1),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(
    grep("Invalid 'phi': must be numeric vector with length 'nrow[(]u0[)]'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(-1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(
    grep("Invalid 'phi': must be numeric vector with non-negative values",
         res[[1]]$message)) > 0)

## Check missing upsilon
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'upsilon' is missing",
                      res[[1]]$message)) > 0)

## Check missing gamma
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'gamma' is missing",
                      res[[1]]$message)) > 0)

## Check missing alpha
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'alpha' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_t1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t1' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_t2
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t2' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_t3
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t3' is missing",
                      res[[1]]$message)) > 0)

## Check missing beta_t4
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t4' is missing",
                      res[[1]]$message)) > 0)

## Check missing end_t1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' is missing",
                      res[[1]]$message)) > 0)

## Check missing end_t2
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t2' is missing",
                      res[[1]]$message)) > 0)

## Check missing end_t3
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t3' is missing",
                      res[[1]]$message)) > 0)

## Check missing end_t4
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' is missing",
                      res[[1]]$message)) > 0)

## Check missing epsilon
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365))
stopifnot(length(grep("'epsilon' is missing",
                      res[[1]]$message)) > 0)

## Check non-numeric upsilon
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = "0.0357",
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'upsilon' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric gamma
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = "0.1",
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'gamma' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric alpha
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = "1.0",
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'alpha' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_t1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = "0.19",
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t1' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_t2
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = "0.085",
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t2' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_t3
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = "0.075",
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t3' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric beta_t4
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = "0.185",
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t4' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-integer end_t1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = "91",
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91.5,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' must be integer",
                      res[[1]]$message)) > 0)

## Check non-integer end_t2
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = "182",
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t2' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182.5,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t2' must be integer",
                      res[[1]]$message)) > 0)

## Check non-integer end_t3
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = "273",
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t3' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273.5,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t3' must be integer",
                      res[[1]]$message)) > 0)

## Check non-integer end_t4
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = "365",
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' must be integer",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365.5,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' must be integer",
                      res[[1]]$message)) > 0)

## Check non-numeric epsilon
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = "0.000011"))
stopifnot(length(grep("'epsilon' must be numeric",
                      res[[1]]$message)) > 0)

## Check that length of upsilon equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = c(0.0357, 0.0357),
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'upsilon' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of gamma equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = c(0.1, 0.1),
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'gamma' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of alpha equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = c(1.0, 1.0),
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'alpha' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_t1 equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = c(0.19, 0.19),
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t1' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_t2 equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = c(0.085, 0.085),
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t2' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_t3 equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = c(0.075, 0.075),
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t3' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of beta_t4 equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = c(0.185, 0.185),
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'beta_t4' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of end_t1 equals 1 or nrow(u0)
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = c(91, 91),
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' must be of length 1 or 'nrow\\(u0\\)'",
                      res[[1]]$message)) > 0)

## Check that length of end_t2 equals 1 or nrow(u0)
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = c(182, 182),
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t2' must be of length 1 or 'nrow\\(u0\\)'",
                      res[[1]]$message)) > 0)

## Check that length of end_t3 equals 1 or nrow(u0)
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = c(273, 273),
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t3' must be of length 1 or 'nrow\\(u0\\)'",
                      res[[1]]$message)) > 0)

## Check that length of end_t4 equals 1 or nrow(u0)
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = c(365, 365),
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' must be of length 1 or 'nrow\\(u0\\)'",
                      res[[1]]$message)) > 0)

## Check that length of epsilon equals 1
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = c(0.000011, 0.000011)))
stopifnot(length(grep("'epsilon' must be of length 1",
                      res[[1]]$message)) > 0)

## Check interval endpoints
res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = -1,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' must be greater than or equal to '0'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 18,
                               end_t3  = 273,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t1' must be less than 'end_t2'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 173,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t2' must be less than 'end_t3'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 365,
                               end_t4  = 365,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t3' must be less than '364'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = -1,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' must be greater than or equal to '0'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 91,
                               end_t2  = 182,
                               end_t3  = 273,
                               end_t4  = 366,
                               epsilon = 0.000011))
stopifnot(length(grep("'end_t4' must be less than or equal to '365'",
                      res[[1]]$message)) > 0)

res <- tools::assertError(SISe(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               phi     = rep(1, nrow(u0)),
                               upsilon = 0.0357,
                               gamma   = 0.1,
                               alpha   = 1.0,
                               beta_t1 = 0.19,
                               beta_t2 = 0.085,
                               beta_t3 = 0.075,
                               beta_t4 = 0.185,
                               end_t1  = 4:9,
                               end_t2  = 5:10,
                               end_t3  = c(8:12, 16),
                               end_t4  = c(2, 11:15),
                               epsilon = 0.000011))
stopifnot(length(grep(
    "'end_t4' must be less than 'end_t1' or greater than 'end_t3'",
    res[[1]]$message)) > 0)

## Check extraction of data from 'suscpetible', and 'infected'
## compartments
model <- SISe(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              phi     = seq(0, by = 0.1, length.out = nrow(u0)),
              upsilon = 0,
              gamma   = 0.1,
              alpha   = 1.0,
              beta_t1 = 0,
              beta_t2 = 0,
              beta_t3 = 0,
              beta_t4 = 0,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0)

result <- run(model, threads = 1)

S_expected <- structure(c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L,
                          1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
                          2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L,
                          3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
                          4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                        .Dim = c(6L, 10L),
                        .Dimnames = list(c("S", "S", "S", "S", "S", "S"),
                                         c("0", "1", "2", "3", "4", "5",
                                           "6", "7", "8", "9")))
S_observed <- trajectory(result, compartments = "S", as.is = TRUE)
stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L),
                        .Dimnames = list(c("I", "I", "I", "I", "I", "I"),
                                         c("0", "1", "2", "3", "4", "5",
                                           "6", "7", "8", "9")))
I_observed <- trajectory(result, compartments = "I", as.is = TRUE)
stopifnot(identical(I_observed, I_expected))

## Check output from trajectory method
model@tspan <- c(1,2)
result <- run(model, threads = 1)
res <- tools::assertError(trajectory(result, c("S", "V1"), as.is = TRUE))
stopifnot(length(grep("Select either continuous or discrete compartments",
                      res[[1]]$message)) > 0)

stopifnot(identical(class(trajectory(result, c("S", "V1"))$V1), "numeric"))
stopifnot(identical(class(trajectory(result, c("V1"))$V1), "numeric"))

traj_expected <- structure(list(
    Node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L),
    Time = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L),
    S = c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L)),
    .Names = c("Node", "Time", "S"),
    row.names = c(NA, -12L),
    class = "data.frame")
stopifnot(identical(trajectory(result, c("S", "S")), traj_expected))

traj_expected <- structure(list(
    Node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L),
    Time = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L),
    S = c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
    V1 = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0, 0.1, 0.2, 0.3, 0.4, 0.5)),
    .Names = c("Node", "Time", "S", "V1"),
    row.names = c(NA, -12L),
    class = "data.frame")
stopifnot(identical(trajectory(result, c("S", "S", "V1", "V1"))[, -4], traj_expected[, -4]))
stopifnot(identical(trajectory(result, c("V1", "V1", "S", "S"))[, -4], traj_expected[, -4]))
stopifnot(all(abs(trajectory(result, c("V1", "V1", "S", "S"))[, 4] - traj_expected$V1) < 1e-8))

## Check extracting a subset of nodes
traj_expected <- structure(list(
    Node = c(2L, 5L, 2L, 5L),
    Time = c(1L, 1L, 2L, 2L),
    S = c(1L, 4L, 1L, 4L),
    V1 = c(0.1, 0.4, 0.1, 0.4)),
    .Names = c("Node", "Time", "S", "V1"),
    row.names = c(NA, -4L),
    class = "data.frame")
stopifnot(identical(trajectory(result, c("S", "S", "V1", "V1"), i = c(5, 2))[, -4], traj_expected[, -4]))
stopifnot(identical(trajectory(result, c("V1", "V1", "S", "S"), i = c(5, 2))[, -4], traj_expected[, -4]))
stopifnot(all(abs(trajectory(result, c("V1", "V1", "S", "S"), i = c(5, 2))[, 4] - traj_expected$V1) < 1e-8))
stopifnot(identical(trajectory(result, c("S", "V1"), i = c(5, 2)), traj_expected))

## Check extracting all compartments in U
stopifnot(identical(trajectory(result, c("S", "I"), as.is = TRUE), result@U))

## Check extracting all compartments of U in internal format
stopifnot(identical(trajectory(result, c("S", "I"), as.is = TRUE), result@U))

## Check extracting a subset compartments of V in internal format
traj_observed <- trajectory(result, "V1", i = c(5, 2), as.is = TRUE)
stopifnot(identical(dim(traj_observed), c(2L, 2L)))
stopifnot(all(abs(traj_observed[1, ] - 0.1) < 1e-8))
stopifnot(all(abs(traj_observed[2, ] - 0.4) < 1e-8))

## Check extracting all compartments in a subset of U as a data.frame
traj_expected <- data.frame(Node = c(2L, 5L, 2L, 5L),
                            Time = c(1L, 1L, 2L, 2L),
                            S = c(1L, 4L, 1L, 4L),
                            I = c(0L, 0L, 0L, 0L))
stopifnot(identical(trajectory(result, c("S", "I"), i = c(5, 2)), traj_expected))

## Check extracting all compartments in U as a data.frame
traj_expected <- data.frame(Node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L),
                            Time = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L),
                            S = c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                            I = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L))
traj_observed <- trajectory(result, c("S", "I"))
stopifnot(identical(traj_observed, traj_expected))

## Check extracting all compartments in U and V as a data.frame
traj_expected <- data.frame(Node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L),
                            Time = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L),
                            S = c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                            I = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                            V1 = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0, 0.1, 0.2, 0.3, 0.4, 0.5))
traj_observed <- trajectory(result, c("S", "I", "V1"))
stopifnot(identical(traj_observed[, -5], traj_expected[, -5]))
stopifnot(all(abs(traj_observed[, 5] - traj_expected$V1) < 1e-8))

## Check SISe plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SISe run function fails for misspecified SISe model
res <- tools::assertError(.Call("SISe_run", NULL, NULL, NULL, NULL,
                                PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid model.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(.Call("SISe_run", "SISe", NULL, NULL, NULL,
                                PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid model.",
                      res[[1]]$message)) > 0)

## Check error non-finite v
model <- SISe(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              phi     = rep(1, nrow(u0)),
              upsilon = 0.0357,
              gamma   = 0.1,
              alpha   = 1.0,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0.000011)
model@gdata["beta_t1"] <- Inf
model@gdata["beta_t2"] <- Inf
model@gdata["beta_t3"] <- Inf
model@gdata["beta_t4"] <- Inf
res <- tools::assertError(run(model, threads = 1))
stopifnot(length(grep("The continuous state 'v' is not finite.",
                      res[[1]]$message)) > 0)

## Check negative v
model <- SISe(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              phi     = rep(1, nrow(u0)),
              upsilon = 0.0357,
              gamma   = 0.1,
              alpha   = 1.0,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = -10.000011)
res <- tools::assertError(run(model, threads = 1))
stopifnot(length(grep("The continuous state 'v' is negative.",
                      res[[1]]$message)) > 0)

## Check data
stopifnot(identical(nrow(events_SISe()), 466692L))
stopifnot(identical(nrow(u0_SISe()), 1600L))
