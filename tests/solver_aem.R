## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2020 Stefan Widgren
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
library(tools)
source("util/check.R")

## Specify the number of threads to use.
set_num_threads(1)

## For debugging
sessionInfo()

## Check invalid u0
res <- assertError(SISe(u0 = "u0"))
check_error(res, "Missing columns in u0.")

u0 <- data.frame(S  = c(9, 9, 9, 9, 9, 10),
                 I  = c(1, 1, 1, 1, 1, 0))

## Check missing columns in u0
res <- assertError(SISe(u0 = u0[, "I", drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- assertError(SISe(u0 = u0[, "S", drop = FALSE]))
check_error(res, "Missing columns in u0.")

## Check 'susceptible' and 'infected' compartments
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

set.seed(22)
result <- run(model, solver = "aem")

S_expected <- structure(c(
    9L, 9L, 10L, 9L, 9L, 10L, 9L, 9L, 10L, 9L, 9L, 10L, 9L,
    9L, 10L, 9L, 9L, 10L, 9L, 9L, 10L, 9L, 8L, 10L, 9L, 9L,
    10L, 9L, 8L, 10L, 9L, 8L, 10L, 9L, 8L, 10L, 9L, 8L, 10L,
    10L, 8L, 10L, 9L, 8L, 10L, 10L, 7L, 10L, 10L, 7L, 10L,
    10L, 7L, 10L, 10L, 7L, 10L, 10L, 7L, 10L),
    .Dim = c(6L, 10L))

S_observed <- trajectory(result, compartments = "S", format = "matrix")
stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L,
                          1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 2L, 0L, 1L, 1L,
                          0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L,
                          0L, 2L, 0L, 1L, 2L, 0L, 0L, 3L, 0L, 0L, 3L, 0L, 0L,
                          3L, 0L, 0L, 3L, 0L, 0L, 3L, 0L),
                        .Dim = c(6L, 10L))

I_observed <- trajectory(result, compartments = "I", format = "matrix")
stopifnot(identical(I_observed, I_expected))

## test with events.
u0 <- data.frame(S = c(10, 9),
                 I = c(0, 1))

events <- data.frame(event      = c(3, 3),
                     time       = c(1, 5),
                     node       = c(1, 2),
                     dest       = c(2, 1),
                     n          = c(2, 2),
                     proportion = c(0, 0),
                     select     = c(2, 2),
                     shift      = c(0, 0))

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

set.seed(123)
result <- run(model, solver = "aem")

S_expected <- structure(c(10L, 8L, 8L, 9L, 7L, 10L, 6L, 10L, 6L, 10L, 8L, 6L,
                          7L, 7L, 7L, 7L, 7L, 7L, 7L, 9L),
                        .Dim = c(2L, 10L))

S_observed <- trajectory(result, compartments = "S", format = "matrix")
stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 2L, 0L, 3L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 4L, 3L,
                          3L, 3L, 3L, 3L, 3L, 3L, 1L),
                        .Dim = c(2L, 10L))

I_observed <- trajectory(result, compartments = "I", format = "matrix")
stopifnot(identical(I_observed, I_expected))

## run with AEM using multiple threads
if (SimInf:::have_openmp()) {
    set.seed(123)
    set_num_threads(2)
    result <- run(model, solver = "aem")
    set_num_threads(1)
    result

    stopifnot(identical(
        length(trajectory(result, compartments = "S", format = "matrix")),
        20L))
    stopifnot(identical(
        length(trajectory(result, compartments = "I", format = "matrix")),
        20L))

    p <- prevalence(result, I ~ S + I, format = "matrix")
    stopifnot(identical(dim(p), c(1L, 10L)))

    p <- prevalence(result, I ~ S + I, level = 3, format = "matrix")
    stopifnot(identical(dim(p), c(2L, 10L)))
}

## Check solver argument
assertError(run(model, solver = 1))
assertError(run(model, solver = c("ssa", "aem")))
assertError(run(model, solver = NA_character_))
assertError(run(model, solver = "non-existing-solver"))
