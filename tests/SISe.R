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
res <- tools::assertError(SISe(u0 = "u0"))
check_error(res, "Missing columns in u0.")

u0 <- structure(list(S  = c(0, 1, 2, 3, 4, 5),
                     I  = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S", "I"),
                row.names = c(NA, -6L), class = "data.frame")

## Check missing columns in u0
res <- tools::assertError(SISe(u0 = u0[, "I", drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(SISe(u0 = u0[, "S", drop = FALSE]))
check_error(res, "Missing columns in u0.")

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
check_error(res, "Invalid 'phi': must be numeric vector.")

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
check_error(res, "Invalid 'phi': must be numeric vector with non-negative values.")

res <- SISe(u0      = u0,
            tspan   = seq_len(10) - 1,
            events  = NULL,
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
stopifnot(identical(res@v0,
                    structure(c(0, 0, 0, 0, 0, 0),
                              .Dim = c(1L, 6L),
                              .Dimnames = list("phi", NULL))))

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
check_error(res, "'upsilon' is missing.")

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
check_error(res, "'gamma' is missing.")

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
check_error(res, "'alpha' is missing.")

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
check_error(res, "'beta_t1' is missing.")

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
check_error(res, "'beta_t2' is missing.")

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
check_error(res, "'beta_t3' is missing.")

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
check_error(res, "'beta_t4' is missing.")

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
check_error(res, "'end_t1' is missing.")

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
check_error(res, "'end_t2' is missing.")

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
check_error(res, "'end_t3' is missing.")

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
check_error(res, "'end_t4' is missing.")

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
check_error(res, "'epsilon' is missing.")

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
check_error(res, "'upsilon' must be numeric.")

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
check_error(res, "'gamma' must be numeric.")

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
check_error(res, "'alpha' must be numeric.")

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
check_error(res, "'beta_t1' must be numeric.")

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
check_error(res, "'beta_t2' must be numeric.")

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
check_error(res, "'beta_t3' must be numeric.")

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
check_error(res, "'beta_t4' must be numeric.")

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
check_error(res, "'end_t1' must be integer.")

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
check_error(res, "'end_t1' must be integer.")

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
check_error(res, "'end_t2' must be integer.")

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
check_error(res, "'end_t2' must be integer.")

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
check_error(res, "'end_t3' must be integer.")

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
check_error(res, "'end_t3' must be integer.")

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
check_error(res, "'end_t4' must be integer.")

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
check_error(res, "'end_t4' must be integer.")

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
check_error(res, "'epsilon' must be numeric.")

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
check_error(res, "'upsilon' must be of length 1.")

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
check_error(res, "'gamma' must be of length 1.")

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
check_error(res, "'alpha' must be of length 1.")

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
check_error(res, "'beta_t1' must be of length 1.")

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
check_error(res, "'beta_t2' must be of length 1.")

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
check_error(res, "'beta_t3' must be of length 1.")

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
check_error(res, "'beta_t4' must be of length 1.")

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
check_error(res, "'epsilon' must be of length 1.")

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
check_error(res, "'end_t1' must be greater than or equal to '0'.")

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
check_error(res, "'end_t1' must be less than 'end_t2'.")

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
check_error(res, "'end_t2' must be less than 'end_t3' or 'end_t3' less than 'end_t1'.")

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
check_error(res, "'end_t3' must be less than '364'.")

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
check_error(res, "'end_t4' must be greater than or equal to '0'.")

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
check_error(res, "'end_t4' must be less than or equal to '365'.")

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
check_error(res, "'end_t4' must be less than 'end_t1' or greater than 'end_t3'.")

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

result <- run(model)

S_expected <- structure(c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L,
                          1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
                          2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L,
                          3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
                          4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                        .Dim = c(6L, 10L))
S_observed <- trajectory(result, compartments = "S", as.is = TRUE)
stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L))
I_observed <- trajectory(result, compartments = "I", as.is = TRUE)
stopifnot(identical(I_observed, I_expected))

## Check output from trajectory method
model@tspan <- c(1,2)
result <- run(model)
res <- tools::assertError(trajectory(result, c("S", "phi"), as.is = TRUE))
check_error(res, "Select either continuous or discrete compartments.")

stopifnot(identical(class(trajectory(result, c("S", "phi"))$phi), "numeric"))
stopifnot(identical(class(trajectory(result, c("phi"))$phi), "numeric"))

traj_expected <- structure(list(
    node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L),
    time = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L),
    S = c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L)),
    .Names = c("node", "time", "S"),
    row.names = c(NA, -12L),
    class = "data.frame")
stopifnot(identical(trajectory(result, c("S", "S")), traj_expected))

traj_expected <- structure(list(
    node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L),
    time = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L),
    S = c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
    phi = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0, 0.1, 0.2, 0.3, 0.4, 0.5)),
    .Names = c("node", "time", "S", "phi"),
    row.names = c(NA, -12L),
    class = "data.frame")
stopifnot(identical(trajectory(result, c("S", "S", "phi", "phi"))[, -4], traj_expected[, -4]))
stopifnot(identical(trajectory(result, c("phi", "phi", "S", "S"))[, -4], traj_expected[, -4]))
stopifnot(all(abs(trajectory(result, c("phi", "phi", "S", "S"))[, 4] - traj_expected$phi) < 1e-8))

## Check extracting a subset of nodes
traj_expected <- structure(list(
    node = c(2L, 5L, 2L, 5L),
    time = c(1L, 1L, 2L, 2L),
    S = c(1L, 4L, 1L, 4L),
    phi = c(0.1, 0.4, 0.1, 0.4)),
    .Names = c("node", "time", "S", "phi"),
    row.names = c(NA, -4L),
    class = "data.frame")
stopifnot(identical(trajectory(result, c("S", "S", "phi", "phi"), node = c(5, 2))[, -4], traj_expected[, -4]))
stopifnot(identical(trajectory(result, c("phi", "phi", "S", "S"), node = c(5, 2))[, -4], traj_expected[, -4]))
stopifnot(all(abs(trajectory(result, c("phi", "phi", "S", "S"), node = c(5, 2))[, 4] - traj_expected$phi) < 1e-8))
stopifnot(identical(trajectory(result, c("S", "phi"), node = c(5, 2)), traj_expected))

## Check extracting all compartments in U
stopifnot(identical(trajectory(result, c("S", "I"), as.is = TRUE), result@U))

## Check extracting all compartments of U in internal format
stopifnot(identical(trajectory(result, c("S", "I"), as.is = TRUE), result@U))

## Check extracting a subset compartments of V in internal format
traj_observed <- trajectory(result, "phi", node = c(5, 2), as.is = TRUE)
stopifnot(identical(dim(traj_observed), c(2L, 2L)))
stopifnot(all(abs(traj_observed[1, ] - 0.1) < 1e-8))
stopifnot(all(abs(traj_observed[2, ] - 0.4) < 1e-8))

## Check extracting all compartments in a subset of U as a data.frame
traj_expected <- data.frame(node = c(2L, 5L, 2L, 5L),
                            time = c(1L, 1L, 2L, 2L),
                            S = c(1L, 4L, 1L, 4L),
                            I = c(0L, 0L, 0L, 0L))
stopifnot(identical(trajectory(result, c("S", "I"), node = c(5, 2)), traj_expected))

## Check extracting all compartments in U as a data.frame
traj_expected <- data.frame(node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L),
                            time = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L),
                            S = c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                            I = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L))
traj_observed <- trajectory(result, c("S", "I"))
stopifnot(identical(traj_observed, traj_expected))

## Check extracting all compartments in U and V as a data.frame
traj_expected <- data.frame(node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L),
                            time = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L),
                            S = c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                            I = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                            phi = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0, 0.1, 0.2, 0.3, 0.4, 0.5))
traj_observed <- trajectory(result, c("S", "I", "phi"))
stopifnot(identical(traj_observed[, -5], traj_expected[, -5]))
stopifnot(all(abs(traj_observed[, 5] - traj_expected$phi) < 1e-8))

## Use formula notation
traj_observed <- trajectory(result, ~.)
stopifnot(identical(traj_observed[, -5], traj_expected[, -5]))
stopifnot(all(abs(traj_observed[, 5] - traj_expected$phi) < 1e-8))

## Specify invalid formula
res <- tools::assertError(trajectory(result, ~1))
check_error(res, "Non-existing compartment(s) in model: '1'.")

## Check SISe plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SISe run function fails for misspecified SISe model
res <- tools::assertError(.Call("SISe_run", NULL, NULL, NULL, PACKAGE = "SimInf"))
check_error(res, "Invalid model.")

res <- tools::assertError(.Call("SISe_run", "SISe", NULL, NULL, PACKAGE = "SimInf"))
check_error(res, "Invalid model.")

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
res <- tools::assertError(run(model))
check_error(res, "The continuous state 'v' is not finite.")

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
res <- tools::assertError(run(model))
check_error(res, "The continuous state 'v' is negative.")

## Check data
stopifnot(identical(nrow(events_SISe()), 466692L))
stopifnot(identical(nrow(u0_SISe()), 1600L))
