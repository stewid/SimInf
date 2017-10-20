## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2016  Stefan Widgren
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

library(SimInf)

## For debugging
sessionInfo()

## Check invalid u0
res <- tools::assertError(SEIR(u0 = "u0"))
stopifnot(length(grep("'u0' must be a data.frame",
                      res[[1]]$message)) > 0)

u0 <- structure(list(S  = c(0, 1, 2, 3, 4, 5),
                     E  = c(0, 0, 0, 0, 0, 0),
                     I  = c(0, 0, 0, 0, 0, 0),
                     R  = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S", "E", "I", "R"),
                row.names = c(NA, -6L), class = "data.frame")

## Check missing columns in u0
res <- tools::assertError(SEIR(u0 = u0[, c("E", "I", "R"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SEIR(u0 = u0[, c("S", "I", "R"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SEIR(u0 = u0[, c("S", "E", "R"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SEIR(u0 = u0[, c("S", "E", "I"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)

## Check missing beta
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               epsilon = 0.5,
                               gamma   = 0.5))
stopifnot(length(grep("'beta' is missing",
                      res[[1]]$message)) > 0)

## Check missing epsilon
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               gamma   = 0.5))
stopifnot(length(grep("'epsilon' is missing",
                      res[[1]]$message)) > 0)

## Check missing gamma
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = 0.5))
stopifnot(length(grep("'gamma' is missing",
                      res[[1]]$message)) > 0)

## Check non-numeric beta
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = "0.5",
                               epsilon = 0.3,
                               gamma   = 0.1))
stopifnot(length(grep("'beta' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric epsilon
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = "0.3",
                               gamma   = 0.1))
stopifnot(length(grep("'epsilon' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric gamma
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = 0.3,
                               gamma   = "0.1"))
stopifnot(length(grep("'gamma' must be numeric",
                      res[[1]]$message)) > 0)

## Check that length of beta equals 1
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = c(0.5, 0.5),
                               epsilon = 0.3,
                               gamma   = 0.1))
stopifnot(length(grep("'beta' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of epsilon equals 1
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = c(0.3, 0.3),
                               gamma   = 0.1))
stopifnot(length(grep("'epsilon' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of gamma equals 1
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = 0.3,
                               gamma   = c(0.1, 0.1)))
stopifnot(length(grep("'gamma' must be of length 1",
                      res[[1]]$message)) > 0)

## Check 'suscpetible', 'infected' and 'recovered' methods
model <- SEIR(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              beta    = 0,
              epsilon = 0,
              gamma   = 0)

result <- run(model, threads = 1)

S_expected <- structure(c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L,
                          1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
                          2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L,
                          3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
                          4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                        .Dim = c(6L, 10L), .Dimnames = list(NULL, NULL))

S_observed <- susceptible(result)

stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L), .Dimnames = list(NULL, NULL))

I_observed <- infected(result)

stopifnot(identical(I_observed, I_expected))

R_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L), .Dimnames = list(NULL, NULL))

R_observed <- recovered(result)

stopifnot(identical(R_observed, R_expected))

## Check data
stopifnot(identical(nrow(events_SEIR()), 466692L))
stopifnot(identical(nrow(u0_SEIR()), 1600L))

## Check SEIR plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SEIR events plot with no events
model <- SEIR(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              beta    = 0,
              epsilon = 0,
              gamma   = 0)
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(model@events)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SEIR events plot method
model <- SEIR(u0      = u0_SEIR(),
              tspan   = seq_len(365 * 4),
              events  = events_SEIR(),
              beta    = 0,
              epsilon = 0,
              gamma   = 0)
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(model@events)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SEIR run function fails for a misspecified SEIR model
res <- tools::assertError(.Call("SEIR_run", NULL, NULL, NULL, NULL,
                                PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid model.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(.Call("SEIR_run", "SEIR", NULL, NULL, NULL,
                                PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid model.",
                      res[[1]]$message)) > 0)
