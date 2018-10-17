## SimInf, a framework for stochastic disease spread simulations
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
res <- tools::assertError(SIR(u0 = "u0"))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)

u0 <- structure(list(S  = c(0, 1, 2, 3, 4, 5),
                     I  = c(0, 0, 0, 0, 0, 0),
                     R  = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S", "I", "R"),
                row.names = c(NA, -6L), class = "data.frame")

## Check missing columns in u0
res <- tools::assertError(SIR(u0 = u0[, c("I", "R"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SIR(u0 = u0[, c("S", "R"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SIR(u0 = u0[, c("S", "I"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)

## Check missing beta
res <- tools::assertError(SIR(u0     = u0,
                              tspan  = seq_len(10) - 1,
                              events = NULL,
                              gamma  = 0.5))
stopifnot(length(grep("'beta' is missing",
                      res[[1]]$message)) > 0)

## Check missing gamma
res <- tools::assertError(SIR(u0     = u0,
                              tspan  = seq_len(10) - 1,
                              events = NULL,
                              beta   = 0.5))
stopifnot(length(grep("'gamma' is missing",
                      res[[1]]$message)) > 0)

## Check non-numeric beta
res <- tools::assertError(SIR(u0      = u0,
                              tspan   = seq_len(10) - 1,
                              events  = NULL,
                              beta    = "0.5",
                              gamma   = 0.1))
stopifnot(length(grep("'beta' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric gamma
res <- tools::assertError(SIR(u0      = u0,
                              tspan   = seq_len(10) - 1,
                              events  = NULL,
                              beta    = 0.5,
                              gamma   = "0.1"))
stopifnot(length(grep("'gamma' must be numeric",
                      res[[1]]$message)) > 0)

## Check that length of beta equals 1
res <- tools::assertError(SIR(u0      = u0,
                              tspan   = seq_len(10) - 1,
                              events  = NULL,
                              beta    = c(0.5, 0.5),
                              gamma   = 0.1))
stopifnot(length(grep("'beta' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of gamma equals 1
res <- tools::assertError(SIR(u0      = u0,
                              tspan   = seq_len(10) - 1,
                              events  = NULL,
                              beta    = 0.5,
                              gamma   = c(0.1, 0.1)))
stopifnot(length(grep("'gamma' must be of length 1",
                      res[[1]]$message)) > 0)

## Extract data from the 'suscpetible', 'infected' and 'recovered'
## compartments
model <- SIR(u0     = u0,
             tspan  = seq_len(10) - 1,
             events = NULL,
             beta   = 0,
             gamma  = 0)

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

R_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L),
                        .Dimnames = list(c("R", "R", "R", "R", "R", "R"),
                                         c("0", "1", "2", "3", "4", "5",
                                           "6", "7", "8", "9")))

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
res <- tools::assertError(.Call("SIR_run", NULL, NULL, NULL, PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid model.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(.Call("SIR_run", "SIR", NULL, NULL, PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid model.",
                      res[[1]]$message)) > 0)

## Check events method
res <- tools::assertError(events())
stopifnot(length(grep("Missing 'model' argument",
                      res[[1]]$message)) > 0)

res <- tools::assertError(events(5))
stopifnot(length(grep("'model' argument is not a 'SimInf_model",
                      res[[1]]$message)) > 0)

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
stopifnot(identical("Invalid rate detected (non-finite or < 0.0)",
                    res[[1]]$message))
res <- tools::assertError(run(model, solver = "aem"))
stopifnot(identical("Invalid rate detected (non-finite or < 0.0)",
                    res[[1]]$message))
