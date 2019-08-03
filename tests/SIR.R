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

library(SimInf)
source("util/check.R")

## Specify the number of threads to use.
set_num_threads(1)

## For debugging
sessionInfo()

## Check invalid u0
res <- tools::assertError(SIR(u0 = "u0"))
check_error(res, "Missing columns in u0.")

u0 <- data.frame(S  = c(0, 1, 2, 3, 4, 5),
                 I  = c(0, 0, 0, 0, 0, 0),
                 R  = c(0, 0, 0, 0, 0, 0))

## Check missing columns in u0
res <- tools::assertError(SIR(u0 = u0[, c("I", "R"), drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(SIR(u0 = u0[, c("S", "R"), drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- tools::assertError(SIR(u0 = u0[, c("S", "I"), drop = FALSE]))
check_error(res, "Missing columns in u0.")

## Check missing beta
res <- tools::assertError(SIR(u0     = u0,
                              tspan  = seq_len(10) - 1,
                              events = NULL,
                              gamma  = 0.5))
check_error(res, "'beta' is missing.")

## Check missing gamma
res <- tools::assertError(SIR(u0     = u0,
                              tspan  = seq_len(10) - 1,
                              events = NULL,
                              beta   = 0.5))
check_error(res, "'gamma' is missing.")

## Check non-numeric beta
res <- tools::assertError(SIR(u0      = u0,
                              tspan   = seq_len(10) - 1,
                              events  = NULL,
                              beta    = "0.5",
                              gamma   = 0.1))
check_error(res, "'beta' must be numeric.")

## Check non-numeric gamma
res <- tools::assertError(SIR(u0      = u0,
                              tspan   = seq_len(10) - 1,
                              events  = NULL,
                              beta    = 0.5,
                              gamma   = "0.1"))
check_error(res, "'gamma' must be numeric.")

## Check that length of beta equals 1
res <- tools::assertError(SIR(u0      = u0,
                              tspan   = seq_len(10) - 1,
                              events  = NULL,
                              beta    = c(0.5, 0.5),
                              gamma   = 0.1))
check_error(res, "'beta' must be of length 1.")

## Check that length of gamma equals 1
res <- tools::assertError(SIR(u0      = u0,
                              tspan   = seq_len(10) - 1,
                              events  = NULL,
                              beta    = 0.5,
                              gamma   = c(0.1, 0.1)))
check_error(res, "'gamma' must be of length 1.")

## Extract data from the 'suscpetible', 'infected' and 'recovered'
## compartments
model <- SIR(u0     = u0,
             tspan  = seq_len(10) - 1,
             events = NULL,
             beta   = 0,
             gamma  = 0)

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

R_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L))

R_observed <- trajectory(result, compartments = "R", as.is = TRUE)

stopifnot(identical(R_observed, R_expected))

## Check data
stopifnot(identical(nrow(events_SIR()), 466692L))
stopifnot(identical(nrow(u0_SIR()), 1600L))

## Check SIR plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SIR events plot with no events
model <- SIR(u0     = u0,
             tspan  = seq_len(10) - 1,
             events = NULL,
             beta   = 0,
             gamma  = 0)
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(model@events)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SIR plot with tspan as Date vector
model <- SIR(u0     = u0,
             tspan  = as.Date(seq_len(10), origin = "2016-12-31"),
             events = NULL,
             beta   = 0,
             gamma  = 0)
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(run(model))
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SIR events plot method
model <- SIR(u0     = u0_SIR(),
             tspan  = seq_len(365 * 4),
             events = events_SIR(),
             beta   = 0,
             gamma  = 0)
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(model@events)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SIR run function fails for misspecified SIR model
res <- tools::assertError(.Call(SimInf:::SIR_run, NULL, NULL, NULL))
check_error(res, "Invalid model.")

res <- tools::assertError(.Call(SimInf:::SIR_run, "SIR", NULL, NULL))
check_error(res, "Invalid model.")

## Check events method
res <- tools::assertError(events())
check_error(res, "Missing 'model' argument.")

res <- tools::assertError(events(5))
check_error(res, "'model' argument is not a 'SimInf_model'.")

model <- SIR(u0     = u0_SIR(),
             tspan  = seq_len(365 * 4),
             events = events_SIR(),
             beta   = 0,
             gamma  = 0)
stopifnot(is(events(model), "SimInf_events"))

## Check that initialisation raises an invalid rate error.
model <- SIR(u0 = data.frame(S = rep(99, 2), I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = -0.077)

res <- tools::assertError(run(model, solver = "ssm"))
check_error(res, "Invalid rate detected (non-finite or < 0.0).")

res <- tools::assertError(run(model, solver = "aem"))
check_error(res, "Invalid rate detected (non-finite or < 0.0).")
