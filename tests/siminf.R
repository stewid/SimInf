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

## 1 Node
## 1 Age category
## 2 Disease-states: Susceptible & Infected
##
## The individual start in the susceptible state, with a probability
## of becoming infected.
init <- structure(list(id = 0, S = 1, I = 0),
                  .Names = c("id", "S", "I"),
                  row.names = c(NA, -1L),
                  class = "data.frame")

model <- SISe(init,
              tspan                       = 0:1000,
              events                      = NULL,
              initial_infectious_pressure = 1,
              response                    = 1,
              recover                     = 1,
              alpha                       = 1,
              beta_q1                     = 1,
              beta_q2                     = 1,
              beta_q3                     = 1,
              beta_q4                     = 1,
              epsilon                     = 1)

result <- run(model, verbose = 0, seed = 123L)

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(identical(sum(result@U), 1001L))
stopifnot(any(result@U[1,]))
stopifnot(any(result@U[2,]))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected.
##
## At t = 1, all individuals are moved to node = 0.
init <- structure(list(id      = c(0, 1, 2, 3, 4, 5),
                       S_age_1 = c(0, 1, 2, 3, 4, 5),
                       I_age_1 = c(0, 0, 0, 0, 0, 0),
                       S_age_2 = c(0, 1, 2, 3, 4, 5),
                       I_age_2 = c(0, 0, 0, 0, 0, 0),
                       S_age_3 = c(0, 1, 2, 3, 4, 5),
                       I_age_3 = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id",
                      "S_age_1", "I_age_1",
                      "S_age_2", "I_age_2",
                      "S_age_3", "I_age_3"),
                  row.names = c(NA, -6L),
                  class = "data.frame")

events <- structure(list(event  = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         time   = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         node   = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest   = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n      = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         prop   = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
                    .Names = c("event", "time", "select", "node", "dest", "n", "prop"),
                    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(init,
               tspan                       = 0:10,
               events                      = events,
               initial_infectious_pressure = rep(1, 6),
               response_age_1              = 1,
               response_age_2              = 1,
               response_age_3              = 1,
               recover_age_1               = 1,
               recover_age_2               = 1,
               recover_age_3               = 1,
               alpha                       = 1,
               beta_q1                     = 1,
               beta_q2                     = 1,
               beta_q3                     = 1,
               beta_q4                     = 1,
               epsilon                     = 1)

result <- run(model, verbose = 0, seed = 123L)

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(all(apply(result@U[1:6,], 1, any)))
stopifnot(identical(sum(result@U[1:6,11]), 45L))
stopifnot(identical(sum(result@U[,1]), 45L))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected.
##
## No external events
init <- structure(list(id      = c(0, 1, 2, 3, 4, 5),
                       S_age_1 = c(0, 1, 2, 3, 4, 5),
                       I_age_1 = c(0, 0, 0, 0, 0, 0),
                       S_age_2 = c(0, 1, 2, 3, 4, 5),
                       I_age_2 = c(0, 0, 0, 0, 0, 0),
                       S_age_3 = c(0, 1, 2, 3, 4, 5),
                       I_age_3 = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id",
                      "S_age_1", "I_age_1",
                      "S_age_2", "I_age_2",
                      "S_age_3", "I_age_3"),
                  row.names = c(NA, -6L), class = "data.frame")

model <- SISe3(init,
               tspan                       = 0:10,
               events                      = NULL,
               initial_infectious_pressure = rep(1, 6),
               response_age_1              = 1,
               response_age_2              = 1,
               response_age_3              = 1,
               recover_age_1               = 1,
               recover_age_2               = 1,
               recover_age_3               = 1,
               alpha                       = 1,
               beta_q1                     = 1,
               beta_q2                     = 1,
               beta_q3                     = 1,
               beta_q4                     = 1,
               epsilon                     = 1)

result <- run(model, verbose = 0, seed = 123L)

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(all(result@U[1:6,] == 0))
stopifnot(all(apply(result@U[seq(from=8, to=36, by=2),], 1, any)))
stopifnot(identical(sum(result@U[,11]), 45L))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a zero probability
## of becoming infected.
##
## At t = 1, all individuals are moved to node = 0.
init <- structure(list(id      = c(0, 1, 2, 3, 4, 5),
                       S_age_1 = c(0, 1, 2, 3, 4, 5),
                       I_age_1 = c(0, 0, 0, 0, 0, 0),
                       S_age_2 = c(0, 1, 2, 3, 4, 5),
                       I_age_2 = c(0, 0, 0, 0, 0, 0),
                       S_age_3 = c(0, 1, 2, 3, 4, 5),
                       I_age_3 = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id",
                      "S_age_1", "I_age_1",
                      "S_age_2", "I_age_2",
                      "S_age_3", "I_age_3"),
                  row.names = c(NA, -6L), class = "data.frame")

events <- structure(list(event  = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         time   = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         node   = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest   = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n      = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         prop   = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
                    .Names = c("event", "time", "select", "node", "dest", "n", "prop"),
                    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(init,
               tspan                       = 0:10,
               events                      = events,
               initial_infectious_pressure = rep(0, 6),
               response_age_1              = 0,
               response_age_2              = 0,
               response_age_3              = 0,
               recover_age_1               = 1,
               recover_age_2               = 1,
               recover_age_3               = 1,
               alpha                       = 0,
               beta_q1                     = 1,
               beta_q2                     = 1,
               beta_q3                     = 1,
               beta_q4                     = 1,
               epsilon                     = 0)

result <- run(model, verbose = 0, seed = 123L)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                 0L, 5L, 0L, 5L, 0L, 5L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L,
                 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L,
                 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L,
                 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 15L, 0L, 15L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
               .Dim = c(36L, 11L))

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(identical(result@U, U))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## Zero probability of becoming infected.
##
## No individuals at t = 0
## At t = 1, all individuals enter in susceptible state
init <- structure(list(id      = c(0, 1, 2, 3, 4, 5),
                       S_age_1 = c(0, 0, 0, 0, 0, 0),
                       I_age_1 = c(0, 0, 0, 0, 0, 0),
                       S_age_2 = c(0, 0, 0, 0, 0, 0),
                       I_age_2 = c(0, 0, 0, 0, 0, 0),
                       S_age_3 = c(0, 0, 0, 0, 0, 0),
                       I_age_3 = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id",
                      "S_age_1", "I_age_1",
                      "S_age_2", "I_age_2",
                      "S_age_3", "I_age_3"),
                  row.names = c(NA, -6L), class = "data.frame")

events <- structure(list(event  = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         time   = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         node   = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest   = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n      = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         prop   = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
                    .Names = c("event", "time", "select", "node", "dest", "n", "prop"),
                    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(init,
               tspan                       = 0:10,
               events                      = events,
               initial_infectious_pressure = rep(0, 6),
               response_age_1              = 0,
               response_age_2              = 0,
               response_age_3              = 0,
               recover_age_1               = 1,
               recover_age_2               = 1,
               recover_age_3               = 1,
               alpha                       = 0,
               beta_q1                     = 1,
               beta_q2                     = 1,
               beta_q3                     = 1,
               beta_q4                     = 1,
               epsilon                     = 0)

result <- run(model, verbose = 0, seed = 123L)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L,
                 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L,
                 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L,
                 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L,
                 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L),
               .Dim = c(36L, 11L))

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(identical(result@U, U))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a zero probability
## of becoming infected.
##
## At t = 3, all individuals exit.
init <- structure(list(id      = c(0, 1, 2, 3, 4, 5),
                       S_age_1 = c(0, 1, 2, 3, 4, 5),
                       I_age_1 = c(0, 0, 0, 0, 0, 0),
                       S_age_2 = c(0, 1, 2, 3, 4, 5),
                       I_age_2 = c(0, 0, 0, 0, 0, 0),
                       S_age_3 = c(0, 1, 2, 3, 4, 5),
                       I_age_3 = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id",
                      "S_age_1", "I_age_1",
                      "S_age_2", "I_age_2",
                      "S_age_3", "I_age_3"),
                  row.names = c(NA, -6L), class = "data.frame")

events <- structure(list(event  = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         time   = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         node   = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest   = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n      = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         prop   = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
                    .Names = c("event", "time", "select", "node", "dest", "n", "prop"),
                    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(init,
               tspan                       = 0:10,
               events                      = events,
               initial_infectious_pressure = rep(0, 6),
               response_age_1              = 0,
               response_age_2              = 0,
               response_age_3              = 0,
               recover_age_1               = 1,
               recover_age_2               = 1,
               recover_age_3               = 1,
               alpha                       = 0,
               beta_q1                     = 1,
               beta_q2                     = 1,
               beta_q3                     = 1,
               beta_q4                     = 1,
               epsilon                     = 0)

result <- run(model, verbose = 0, seed = 123L)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
               .Dim = c(36L, 11L))

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(identical(result@U, U))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a zero probability of
## becoming infected.
##
## At t = 3, all individuals in age_1 category age.
## At t = 6, all individuals in age_2 category age.
init <- structure(list(id      = c(0, 1, 2, 3, 4, 5),
                       S_age_1 = c(0, 1, 2, 3, 4, 5),
                       I_age_1 = c(0, 0, 0, 0, 0, 0),
                       S_age_2 = c(0, 1, 2, 3, 4, 5),
                       I_age_2 = c(0, 0, 0, 0, 0, 0),
                       S_age_3 = c(0, 1, 2, 3, 4, 5),
                       I_age_3 = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id",
                      "S_age_1", "I_age_1",
                      "S_age_2", "I_age_2",
                      "S_age_3", "I_age_3"),
                  row.names = c(NA, -6L), class = "data.frame")

events <- structure(list(event  = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
                         time   = c(3, 3, 3, 3, 3, 6, 6, 6, 6, 6),
                         select = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
                         node   = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
                         dest   = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n      = c(1, 2, 3, 4, 5, 2, 4, 6, 8, 10),
                         prop   = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
                    .Names = c("event", "time", "select", "node", "dest", "n", "prop"),
                    row.names = c(NA, -10L), class = "data.frame")

model <- SISe3(init,
               tspan                       = 0:10,
               events                      = events,
               initial_infectious_pressure = rep(0, 6),
               response_age_1              = 0,
               response_age_2              = 0,
               response_age_3              = 0,
               recover_age_1               = 1,
               recover_age_2               = 1,
               recover_age_3               = 1,
               alpha                       = 0,
               beta_q1                     = 1,
               beta_q2                     = 1,
               beta_q3                     = 1,
               beta_q4                     = 1,
               epsilon                     = 0)

result <- run(model, verbose = 0, seed = 123L)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L,
                 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L,
                 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 1L,
                 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 0L, 3L, 0L, 4L,
                 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 1L, 0L, 1L, 0L, 2L, 0L, 2L, 0L, 2L, 0L, 3L, 0L, 3L,
                 0L, 3L, 0L, 4L, 0L, 4L, 0L, 4L, 0L, 5L, 0L, 5L, 0L, 5L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 1L, 0L, 0L, 0L, 4L, 0L, 2L,
                 0L, 0L, 0L, 6L, 0L, 3L, 0L, 0L, 0L, 8L, 0L, 4L, 0L, 0L, 0L, 10L,
                 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 1L, 0L, 0L,
                 0L, 4L, 0L, 2L, 0L, 0L, 0L, 6L, 0L, 3L, 0L, 0L, 0L, 8L, 0L, 4L,
                 0L, 0L, 0L, 10L, 0L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 2L, 0L, 1L, 0L, 0L, 0L, 4L, 0L, 2L, 0L, 0L, 0L, 6L, 0L, 3L, 0L,
                 0L, 0L, 8L, 0L, 4L, 0L, 0L, 0L, 10L, 0L, 5L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 6L, 0L, 0L,
                 0L, 0L, 0L, 9L, 0L, 0L, 0L, 0L, 0L, 12L, 0L, 0L, 0L, 0L, 0L,
                 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L,
                 0L, 0L, 0L, 6L, 0L, 0L, 0L, 0L, 0L, 9L, 0L, 0L, 0L, 0L, 0L, 12L,
                 0L, 0L, 0L, 0L, 0L, 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 6L, 0L, 0L, 0L, 0L, 0L, 9L, 0L,
                 0L, 0L, 0L, 0L, 12L, 0L, 0L, 0L, 0L, 0L, 15L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L, 0L, 0L, 0L, 6L, 0L, 0L,
                 0L, 0L, 0L, 9L, 0L, 0L, 0L, 0L, 0L, 12L, 0L, 0L, 0L, 0L, 0L,
                 15L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3L, 0L, 0L,
                 0L, 0L, 0L, 6L, 0L, 0L, 0L, 0L, 0L, 9L, 0L, 0L, 0L, 0L, 0L, 12L,
                 0L, 0L, 0L, 0L, 0L, 15L, 0L),
               .Dim = c(36L, 11L))

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(identical(result@U, U))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

## 6 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## All individuals start in susceptible state, with a probability of
## becoming infected and then return to susceptible.
##
## No external events
init <- structure(list(id      = c(0, 1, 2, 3, 4, 5),
                       S_age_1 = c(0, 1, 2, 3, 4, 5),
                       I_age_1 = c(0, 0, 0, 0, 0, 0),
                       S_age_2 = c(0, 1, 2, 3, 4, 5),
                       I_age_2 = c(0, 0, 0, 0, 0, 0),
                       S_age_3 = c(0, 1, 2, 3, 4, 5),
                       I_age_3 = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id",
                      "S_age_1", "I_age_1",
                      "S_age_2", "I_age_2",
                      "S_age_3", "I_age_3"),
                  row.names = c(NA, -6L), class = "data.frame")

model <- SISe3(init,
               tspan                       = 0:10,
               events                      = NULL,
               initial_infectious_pressure = rep(1, 6),
               response_age_1              = 1,
               response_age_2              = 1,
               response_age_3              = 1,
               recover_age_1               = 1,
               recover_age_2               = 1,
               recover_age_3               = 1,
               alpha                       = 1,
               beta_q1                     = 1,
               beta_q2                     = 1,
               beta_q3                     = 1,
               beta_q4                     = 1,
               epsilon                     = 1)

result <- run(model, verbose = 0, seed = 123L)

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(identical(sum(result@U[1:6,]), 0L))
stopifnot(all(apply(result@U[7:36,], 1, any)))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))
