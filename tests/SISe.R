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

init <- structure(list(id = c(0, 1, 2, 3, 4, 5),
                       S  = c(0, 1, 2, 3, 4, 5),
                       I  = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id", "S", "I"),
                  row.names = c(NA, -6L), class = "data.frame")

## Check missing columns in init
tools::assertError(SISe(init = init[, c("S", "I")]))
tools::assertError(SISe(init = init[, c("id", "I")]))
tools::assertError(SISe(init = init[, c("id", "S")]))

## Check missing upsilon
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check missing gamma
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check missing alpha
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check missing beta_q1
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check missing beta_q2
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check missing beta_q3
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check missing beta_q4
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        epsilon = 0.000011))

## Check missing epsilon
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185))

## Check non-numeric upsilon
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = "0.0357",
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check non-numeric gamma
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = "0.1",
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check non-numeric alpha
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = "1.0",
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check non-numeric beta_q1
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = "0.19",
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check non-numeric beta_q2
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = "0.085",
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check non-numeric beta_q3
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = "0.075",
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check non-numeric beta_q4
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = "0.185",
                        epsilon = 0.000011))

## Check non-numeric epsilon
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = "0.000011"))

## Check that length of upsilon equals 1
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = c(0.0357, 0.0357),
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check that length of gamma equals 1
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = c(0.1, 0.1),
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check that length of alpha equals 1
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = c(1.0, 1.0),
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check that length of epsilon equals 1
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = c(0.000011, 0.000011)))

## Check that length of beta_q1 equals 1 or nrow(init)
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = rep(0.19, nrow(init) + 1),
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check that length of beta_q2 equals 1 or nrow(init)
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = rep(0.085, nrow(init) + 1),
                        beta_q3 = 0.075,
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check that length of beta_q3 equals 1 or nrow(init)
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = rep(0.075, nrow(init) + 1),
                        beta_q4 = 0.185,
                        epsilon = 0.000011))

## Check that length of beta_q4 equals 1 or nrow(init)
tools::assertError(SISe(init    = init,
                        tspan   = seq_len(10) - 1,
                        events  = NULL,
                        phi     = rep(1, nrow(init)),
                        upsilon = 0.0357,
                        gamma   = 0.1,
                        alpha   = 1.0,
                        beta_q1 = 0.19,
                        beta_q2 = 0.085,
                        beta_q3 = 0.075,
                        beta_q4 = rep(0.185, nrow(init) + 1),
                        epsilon = 0.000011))

## Check 'suscpetible' and 'infected' methods
model <- SISe(init    = init,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              phi     = rep(0, nrow(init)),
              upsilon = 0.0357,
              gamma   = 0.1,
              alpha   = 1.0,
              beta_q1 = 0.19,
              beta_q2 = 0.085,
              beta_q3 = 0.075,
              beta_q4 = 0.185,
              epsilon = 0.000011)

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

## Check SISe plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SISe run function fails for misspecified SISe model
tools::assertError(.Call(siminf:::SISe_run, NULL, NULL, NULL))
tools::assertError(.Call(siminf:::SISe_run, "SISe", NULL, NULL))

setClass("DummySISe", slots = c(a = "character"))
model <- new("DummySISe", a = "SISe")
tools::assertError(.Call(siminf:::SISe_run, model, NULL, NULL))
