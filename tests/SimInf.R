## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2021 Stefan Widgren
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
max_threads <- set_num_threads(1)

## For debugging
sessionInfo()

## 1 Node
## 1 Age category
## 2 Disease-states: Susceptible & Infected
##
## The individual start in the susceptible state, with a probability
## of becoming infected.
u0 <- data.frame(S = 1, I = 0)

model <- SISe(u0      = u0,
              tspan   = 0:1000,
              events  = NULL,
              phi     = 1,
              upsilon = 1,
              gamma   = 1,
              alpha   = 1,
              beta_t1 = 1,
              beta_t2 = 1,
              beta_t3 = 1,
              beta_t4 = 1,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 1)

result <- run(model)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(sum(trajectory(result, format = "matrix")), 1001L))
stopifnot(any(trajectory(result, format = "matrix")[1, ]))
stopifnot(any(trajectory(result, format = "matrix")[2, ]))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp() && max_threads > 1) {
    set_num_threads(2)
    result_omp <- run(model)
    set_num_threads(1)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(sum(result_omp@U), 1001L))
    stopifnot(any(result_omp@U[1, ]))
    stopifnot(any(result_omp@U[2, ]))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected.
##
## At t = 1, all individuals are moved to node = 1.
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

events <- data.frame(
    event      = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node       = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6),
    shift      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(1, 6),
               upsilon_1 = 1,
               upsilon_2 = 1,
               upsilon_3 = 1,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 1,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 1)

set.seed(123)
result <- run(model)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(all(apply(trajectory(result, format = "matrix")[1:6, ], 1, any)))
stopifnot(identical(sum(trajectory(result, format = "matrix")[1:6, 11]), 45L))
stopifnot(identical(sum(trajectory(result, format = "matrix")[, 1]), 45L))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp() && max_threads > 1) {
    set.seed(123)
    set_num_threads(2)
    result_omp <- run(model)
    set_num_threads(1)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(all(apply(result_omp@U[1:6, ], 1, any)))
    stopifnot(identical(sum(result_omp@U[1:6, 11]), 45L))
    stopifnot(identical(sum(result_omp@U[, 1]), 45L))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected.
##
## At t = 1, all individuals in age category 1 and 2 are moved to node
## = 1 and age.
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

events <- data.frame(
    event      = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node       = c(2, 2, 3, 3, 4, 4, 5, 5, 6, 6),
    dest       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n          = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(4, 5, 4, 5, 4, 5, 4, 5, 4, 5),
    shift      = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2))

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(1, 6),
               upsilon_1 = 1,
               upsilon_2 = 1,
               upsilon_3 = 1,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 1,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 1)

result <- run(model)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))

m <- trajectory(result, compartments = "S_2", index = 1, format = "matrix") +
    trajectory(result, compartments = "I_2", index = 1, format = "matrix")
dimnames(m) <- NULL
stopifnot(identical(m,
                    structure(c(0L, 15L, 15L, 15L, 15L,
                                15L, 15L, 15L, 15L, 15L, 15L),
                              .Dim = c(1L, 11L))))
m <- trajectory(result, compartments = "S_3", index = 1, format = "matrix") +
    trajectory(result, compartments = "I_3", index = 1, format = "matrix")
dimnames(m) <- NULL
stopifnot(identical(m,
                    structure(c(0L, 15L, 15L, 15L, 15L,
                                15L, 15L, 15L, 15L, 15L, 15L),
                              .Dim = c(1L, 11L))))
stopifnot(identical(sum(trajectory(result, format = "matrix")[, 1]), 45L))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp() && max_threads > 1) {
    set_num_threads(2)
    result_omp <- run(model)
    set_num_threads(1)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    m <- trajectory(result,
                    compartments = "S_2",
                    index = 1,
                    format = "matrix") +
        trajectory(result,
                   compartments = "I_2",
                   index = 1,
                   format = "matrix")
    dimnames(m) <- NULL
    stopifnot(identical(m,
                        structure(c(0L, 15L, 15L, 15L, 15L,
                                    15L, 15L, 15L, 15L, 15L, 15L),
                                  .Dim = c(1L, 11L))))
    m <- trajectory(result,
                    compartments = "S_3",
                    index = 1,
                    format = "matrix") +
        trajectory(result,
                   compartments = "I_3",
                   index = 1,
                   format = "matrix")
    dimnames(m) <- NULL
    stopifnot(identical(m,
                        structure(c(0L, 15L, 15L, 15L, 15L,
                                    15L, 15L, 15L, 15L, 15L, 15L),
                                  .Dim = c(1L, 11L))))
    stopifnot(identical(sum(result_omp@U[, 1]), 45L))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected.
##
## No scheduled events
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = NULL,
               phi       = rep(1, 6),
               upsilon_1 = 1,
               upsilon_2 = 1,
               upsilon_3 = 1,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 1,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 1)

result <- run(model)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(all(trajectory(result, format = "matrix")[1:6, ] == 0))
stopifnot(identical(sum(trajectory(result, format = "matrix")[, 11]), 45L))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

i <- seq(from = 8, to = 36, by = 2)
stopifnot(all(apply(trajectory(result, format = "matrix")[i, ], 1, any)))

if (SimInf:::have_openmp() && max_threads > 1) {
    set_num_threads(2)
    result_omp <- run(model)
    set_num_threads(1)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(all(result_omp@U[1:6, ] == 0))
    stopifnot(identical(sum(result_omp@U[, 11]), 45L))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
    stopifnot(all(apply(result_omp@U[i, ], 1, any)))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a zero probability
## of becoming infected.
##
## At t = 1, all individuals are moved to node = 1.
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

events <- data.frame(
    event      = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node       = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6),
    shift      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 6),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

U_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L,
                          0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L,
                          3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L,
                          0L, 5L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L,
                          0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L,
                          15L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L, 0L,
                          15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L,
                          0L, 15L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L,
                          0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(36L, 11L))

result <- run(model)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(trajectory(result, format = "matrix"), U_expected))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp() && max_threads > 1) {
    set_num_threads(2)
    result_omp <- run(model)
    set_num_threads(1)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(result_omp@U, U_expected))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## Zero probability of becoming infected.
##
## No individuals at t = 0
## At t = 1, all individuals enter in susceptible state
u0 <- data.frame(S_1 = c(0, 0, 0, 0, 0, 0),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 0, 0, 0, 0, 0),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 0, 0, 0, 0, 0),
                 I_3 = c(0, 0, 0, 0, 0, 0))

events <- data.frame(
    event      = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node       = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3),
    shift      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 6),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

U_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L,
                          1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L,
                          0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L,
                          5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L,
                          2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L,
                          0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                          0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L,
                          4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L,
                          1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                          0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L,
                          5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L,
                          0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L,
                          3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                          0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L,
                          0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L,
                          4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L,
                          2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L,
                          0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L,
                          5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                          0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L,
                          3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L,
                          0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L,
                          0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L,
                          4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L,
                          2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                          0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L),
                        .Dim = c(36L, 11L))

result <- run(model)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(trajectory(result, format = "matrix"), U_expected))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp() && max_threads > 1) {
    set_num_threads(2)
    result_omp <- run(model)
    set_num_threads(1)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(result_omp@U, U_expected))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a zero probability
## of becoming infected.
##
## At t = 3, all individuals exit.
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

events <- data.frame(
    event      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    time       = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    node       = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6),
    shift      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 6),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

U_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L,
                          0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L,
                          3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L,
                          0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L,
                          1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L,
                          0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L,
                          5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L,
                          2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L,
                          0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(36L, 11L))

result <- run(model)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(trajectory(result, format = "matrix"), U_expected))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp() && max_threads > 1) {
    set_num_threads(2)
    result_omp <- run(model)
    set_num_threads(1)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(result_omp@U, U_expected))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a zero probability of
## becoming infected.
##
## At t = 3, all individuals in age category 1 age.
## At t = 6, all individuals in age category 2 age.
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

events <- data.frame(event      = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
                     time       = c(3, 3, 3, 3, 3, 6, 6, 6, 6, 6),
                     node       = c(2, 3, 4, 5, 6, 2, 3, 4, 5, 6),
                     dest       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                     n          = c(1, 2, 3, 4, 5, 2, 4, 6, 8, 10),
                     proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                     select     = c(4, 4, 4, 4, 4, 5, 5, 5, 5, 5),
                     shift      = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2))

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 6),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0)

U_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L,
                          0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L,
                          3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L,
                          0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L,
                          1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L,
                          0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L,
                          5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L,
                          2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L,
                          0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 1L, 0L, 0L,
                          0L, 4L, 0L, 2L, 0L, 0L, 0L, 6L, 0L, 3L, 0L,
                          0L, 0L, 8L, 0L, 4L, 0L, 0L, 0L, 10L, 0L, 5L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L,
                          1L, 0L, 0L, 0L, 4L, 0L, 2L, 0L, 0L, 0L, 6L,
                          0L, 3L, 0L, 0L, 0L, 8L, 0L, 4L, 0L, 0L, 0L,
                          10L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 2L, 0L, 1L, 0L, 0L, 0L, 4L, 0L, 2L, 0L,
                          0L, 0L, 6L, 0L, 3L, 0L, 0L, 0L, 8L, 0L, 4L,
                          0L, 0L, 0L, 10L, 0L, 5L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L, 0L, 0L,
                          0L, 6L, 0L, 0L, 0L, 0L, 0L, 9L, 0L, 0L, 0L,
                          0L, 0L, 12L, 0L, 0L, 0L, 0L, 0L, 15L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L,
                          0L, 0L, 0L, 0L, 0L, 6L, 0L, 0L, 0L, 0L, 0L,
                          9L, 0L, 0L, 0L, 0L, 0L, 12L, 0L, 0L, 0L, 0L,
                          0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 6L, 0L, 0L,
                          0L, 0L, 0L, 9L, 0L, 0L, 0L, 0L, 0L, 12L, 0L,
                          0L, 0L, 0L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L,
                          6L, 0L, 0L, 0L, 0L, 0L, 9L, 0L, 0L, 0L, 0L,
                          0L, 12L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L,
                          0L, 0L, 0L, 0L, 6L, 0L, 0L, 0L, 0L, 0L, 9L,
                          0L, 0L, 0L, 0L, 0L, 12L, 0L, 0L, 0L, 0L, 0L,
                          15L, 0L),
                        .Dim = c(36L, 11L))

result <- run(model)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(trajectory(result, format = "matrix"), U_expected))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp() && max_threads > 1) {
    set_num_threads(2)
    result_omp <- run(model)
    set_num_threads(1)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(result_omp@U, U_expected))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected and then return to susceptible.
##
## No scheduled events
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = NULL,
               phi       = rep(1, 6),
               upsilon_1 = 1,
               upsilon_2 = 1,
               upsilon_3 = 1,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 1,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 1)

set.seed(123)
result <- run(model)
stopifnot(identical(model@G, result@G))
stopifnot(identical(model@S, result@S))
stopifnot(identical(sum(trajectory(result, format = "matrix")[1:6, ]), 0L))
stopifnot(all(apply(trajectory(result, format = "matrix")[7:36, ], 1, any)))
stopifnot(identical(model@ldata, result@ldata))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

if (SimInf:::have_openmp() && max_threads > 1) {
    set.seed(123)
    set_num_threads(2)
    result_omp <- run(model)
    set_num_threads(1)
    stopifnot(identical(model@G, result_omp@G))
    stopifnot(identical(model@S, result_omp@S))
    stopifnot(identical(sum(result_omp@U[1:6, ]), 0L))
    stopifnot(all(apply(result_omp@U[7:36, ], 1, any)))
    stopifnot(identical(model@ldata, result_omp@ldata))
    stopifnot(identical(model@tspan, result_omp@tspan))
    stopifnot(identical(model@u0, result_omp@u0))
    stopifnot(identical(model@events, result_omp@events))
}

## Check extraction of number of threads
model <- SISe(u0      = data.frame(S = 99, I = 1),
              tspan   = seq_len(1000) - 1,
              events  = NULL,
              phi     = 1,
              upsilon = 0.017,
              gamma   = 0.1,
              alpha   = 1,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0.000011)

.Call(SimInf:::SISe_run, model, NULL)

res <- assertError(set_num_threads(-1L))
check_error(res, "'threads' must be an integer >= 1.")

res <- assertError(set_num_threads(-1))
check_error(res, "'threads' must be an integer >= 1.")

res <- assertError(set_num_threads("1"))
check_error(res, "'threads' must be an integer >= 1.")

res <- assertError(set_num_threads(c(1L, 1L)))
check_error(res, "'threads' must be an integer >= 1.")

res <- assertError(set_num_threads(c(1, 1)))
check_error(res, "'threads' must be an integer >= 1.")

res <- assertError(set_num_threads(NA_integer_))
check_error(res, "'threads' must be an integer >= 1.")

res <- assertError(set_num_threads(NA_real_))
check_error(res, "'threads' must be an integer >= 1.")

## Check that we have at least one available thread.  First, try to
## set the number of threads to 0 and then set the number of threads
## to one to get the previous number of threads.
.Call(SimInf:::SimInf_init_threads, 0L)
stopifnot(identical(.Call(SimInf:::SimInf_init_threads, 1L), 1L))
