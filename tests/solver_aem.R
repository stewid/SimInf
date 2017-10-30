## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2017  Robin Eriksson
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

library(SimInf)

## For debugging
sessionInfo()

## Check invalid u0
res <- tools::assertError(SISe(u0 = "u0"))
stopifnot(length(grep("'u0' must be a data.frame",
                      res[[1]]$message)) > 0)

u0 <- structure(list(S  = c(9, 9, 9, 9, 9, 10),
                     I  = c(1, 1, 1, 1, 1, 0)),
                .Names = c("S", "I"),
                row.names = c(NA, -6L), class = "data.frame")

## Check missing columns in u0
res <- tools::assertError(SISe(u0 = u0[, "I", drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SISe(u0 = u0[, "S", drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)

## Check 'susceptible' and 'infected' methods
## no events
model <- SISe(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              phi     = rep(0, nrow(u0)),
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

result <- run(model, seed = 22, threads = 1, solver = "aem")


S_expected <- structure(c(9L, 9L, 9L, 9L, 9L, 10L, 9L, 10L, 9L, 9L, 10L, 10L,
                          9L, 10L, 9L, 9L, 10L, 10L, 9L, 10L, 9L, 9L, 10L, 10L,
                          9L, 10L, 9L, 9L, 10L, 10L, 9L, 10L, 9L, 9L, 10L, 10L,
                          9L, 10L, 8L, 9L, 10L, 10L, 9L, 10L, 7L, 9L, 10L, 10L,
                          9L, 10L, 6L, 9L, 10L, 10L, 8L, 10L, 6L, 8L, 10L, 10L),
                        .Dim = c(6L, 10L), .Dimnames = list(NULL, NULL))

S_observed <- susceptible(result)

stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(1L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 1L,
                          0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L,
                          1L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 2L,
                          1L, 0L, 0L, 1L, 0L, 3L, 1L, 0L, 0L, 1L, 0L, 4L, 1L,
                          0L, 0L, 2L, 0L, 4L, 2L, 0L, 0L),
                        .Dim = c(6L, 10L), .Dimnames = list(NULL, NULL))

I_observed <- infected(result)

stopifnot(identical(I_observed, I_expected))

## test with events.
u0 <- structure(list(S = c(10, 9),
                     I = c(0, 1)),
                .Names = c("S", "I"),
                row.names = c(NA, -2L),
                class = "data.frame")

events <- structure(list(
    event      = c(3, 3),
    time       = c(1, 5),
    node       = c(1, 2),
    dest       = c(2, 1),
    n          = c(2, 2),
    proportion = c(0, 0),
    select     = c(2, 2),
    shift      = c(0, 0)),
    .Names = c("event", "time", "node", "dest",
               "n", "proportion", "select", "shift"),
    row.names = c(NA, -2L), class = "data.frame")


model <- SISe(u0  = u0,
              tspan   = seq_len(10) - 1,
              events  = events,
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

result <- run(model, threads = 1, seed = 123L, solver = "aem")


S_expected <- structure(c(10L, 9L, 8L, 9L, 8L, 9L, 7L, 8L, 7L, 8L, 10L,
                          6L, 10L, 6L, 10L, 6L, 10L, 5L, 10L, 5L),
                        .Dim = c(2L, 10L), .Dimnames = list(NULL, NULL))

S_observed <- susceptible(result)

stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 1L, 0L, 3L, 0L, 3L, 1L, 4L, 1L, 4L, 0L, 4L,
                          0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L),
                        .Dim = c(2L, 10L), .Dimnames = list(NULL, NULL))

I_observed <- infected(result)

stopifnot(identical(I_observed, I_expected))

## run with AEM using multiple threads
if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 123L, solver = "aem")
    result_omp

    stopifnot(identical(length(susceptible(result_omp)), 20L))
    stopifnot(identical(length(infected(result_omp)), 20L))
    stopifnot(identical(length(prevalence(result_omp)), 10L))
    stopifnot(is.null(dim(prevalence(result_omp))))
    stopifnot(identical(dim(prevalence(result_omp, type = "wnp")), c(2L, 10L)))
}

## Check solver argument
res <- tools::assertError(run(model, threads = 1, solver = 1))
stopifnot(length(grep("Invalid 'solver' value.",
                      res[[1]]$message)) > 0)
res <- tools::assertError(run(model, threads = 1, solver = c("ssa", "aem")))
stopifnot(length(grep("Invalid 'solver' value.",
                      res[[1]]$message)) > 0)
res <- tools::assertError(run(model, threads = 1, solver = NA_character_))
stopifnot(length(grep("Invalid 'solver' value.",
                      res[[1]]$message)) > 0)
res <- tools::assertError(run(model, threads = 1, solver = "non-existing-solver"))
stopifnot(length(grep("Invalid 'solver' value.",
                      res[[1]]$message)) > 0)
