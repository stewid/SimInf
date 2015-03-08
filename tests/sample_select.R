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

## Test that sample_select in events.c works

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 1, with a zero
## probability of becoming infected.
##
## At t = 1, two individuals are moved to node = 0. This should fail.
init <- structure(list(id  = c(0, 1),
                       S_1 = c(0, 1),
                       I_1 = c(0, 0),
                       S_2 = c(0, 0),
                       I_2 = c(0, 0),
                       S_3 = c(0, 0),
                       I_3 = c(0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -2L),
                  class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         select     = 0,
                         node       = 1,
                         dest       = 0,
                         n          = 2,
                         proportion = 1),
                    .Names = c("event", "time", "select", "node", "dest", "n", "proportion"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(init,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_q1   = 1,
               beta_q2   = 1,
               beta_q3   = 1,
               beta_q4   = 1,
               epsilon   = 0)

tools::assertError(run(model, verbose = 0))

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 1, with a zero
## probability of becoming infected.
##
## At t = 1, -1 individuals are moved to node = 0. This should fail.
init <- structure(list(id  = c(0, 1),
                       S_1 = c(0, 1),
                       I_1 = c(0, 0),
                       S_2 = c(0, 0),
                       I_2 = c(0, 0),
                       S_3 = c(0, 0),
                       I_3 = c(0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         select     = 0,
                         node       = 1,
                         dest       = 0,
                         n          = -1,
                         proportion = 1),
                    .Names = c("event", "time", "select", "node", "dest", "n", "proportion"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(init,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_q1   = 1,
               beta_q2   = 1,
               beta_q3   = 1,
               beta_q4   = 1,
               epsilon   = 0)

tools::assertError(run(model, verbose = 0))

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 1, with a zero
## probability of becoming infected.
##
## At t = 1, a proportion of 10 individuals are moved to node =
## 0. This should fail.
init <- structure(list(id  = c(0, 1),
                       S_1 = c(0, 1),
                       I_1 = c(0, 0),
                       S_2 = c(0, 0),
                       I_2 = c(0, 0),
                       S_3 = c(0, 0),
                       I_3 = c(0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         select     = 0,
                         node       = 1,
                         dest       = 0,
                         n          = 0,
                         proportion = 10),
                    .Names = c("event", "time", "select", "node", "dest", "n", "proportion"),
                    row.names = c(NA, -1L), class = "data.frame")

## We should not be able to create model with prop = 10
tools::assertError(SISe3(init,
                         tspan     = 0:10,
                         events    = events,
                         phi       = rep(0, 2),
                         upsilon_1 = 0,
                         upsilon_2 = 0,
                         upsilon_3 = 0,
                         gamma_1   = 1,
                         gamma_2   = 1,
                         gamma_3   = 1,
                         alpha     = 0,
                         beta_q1   = 1,
                         beta_q2   = 1,
                         beta_q3   = 1,
                         beta_q4   = 1,
                         epsilon   = 0))

## Replace proportion = 10 to proportion = 1
events$proportion <- 1

model <- SISe3(init,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_q1   = 1,
               beta_q2   = 1,
               beta_q3   = 1,
               beta_q4   = 1,
               epsilon   = 0)

## Replace proportion = 10 to proportion = 1
model@events@proportion <- 10

tools::assertError(run(model, verbose = 0))

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 1, with a zero
## probability of becoming infected.
##
## At t = 1, a proportion of -1 individuals are moved to node =
## 0. This should fail.
init <- structure(list(id  = c(0, 1),
                       S_1 = c(0, 1),
                       I_1 = c(0, 0),
                       S_2 = c(0, 0),
                       I_2 = c(0, 0),
                       S_3 = c(0, 0),
                       I_3 = c(0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         select     = 0,
                         node       = 1,
                         dest       = 0,
                         n          = 0,
                         proportion = -1),
                    .Names = c("event", "time", "select", "node", "dest", "n", "proportion"),
                    row.names = c(NA, -1L), class = "data.frame")

## We should not be able to create model with proportion = -1
tools::assertError(SISe3(init,
                         tspan     = 0:10,
                         events    = events,
                         phi       = rep(0, 2),
                         upsilon_1 = 0,
                         upsilon_2 = 0,
                         upsilon_3 = 0,
                         gamma_1   = 1,
                         gamma_2   = 1,
                         gamma_3   = 1,
                         alpha     = 0,
                         beta_q1   = 1,
                         beta_q2   = 1,
                         beta_q3   = 1,
                         beta_q4   = 1,
                         epsilon   = 0))

## Replace proportion = -1 to proportion = 0
events$proportion <- 0

model <- SISe3(init,
               tspan     = 0:10,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_q1   = 1,
               beta_q2   = 1,
               beta_q3   = 1,
               beta_q4   = 1,
               epsilon   = 0)

## Replace proportion = 0 to proportion = -1
model@events@proportion <- -1

tools::assertError(run(model, verbose = 0))

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 1, with a zero
## probability of becoming infected.
##
## At t = 1, a proportion of 0 individuals are moved to node = 0.
init <- structure(list(id  = c(0, 1),
                       S_1 = c(0, 1),
                       I_1 = c(0, 0),
                       S_2 = c(0, 0),
                       I_2 = c(0, 0),
                       S_3 = c(0, 0),
                       I_3  = c(0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         select     = 0,
                         node       = 1,
                         dest       = 0,
                         n          = 0,
                         proportion = 0),
                    .Names = c("event", "time", "select", "node", "dest", "n", "proportion"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(init,
               tspan     = 0:2,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_q1   = 1,
               beta_q2   = 1,
               beta_q3   = 1,
               beta_q4   = 1,
               epsilon   = 0)

result <- run(model, verbose = 0)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 0L, 0L, 0L, 0L), .Dim = c(12L, 3L))

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(identical(result@U, U))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## One individual start in susceptible state in node = 1, with a zero
## probability of becoming infected.
##
## At t = 1, proportion of all (1) individuals are moved to node = 0.
init <- structure(list(id  = c(0, 1),
                       S_1 = c(0, 1),
                       I_1 = c(0, 0),
                       S_2 = c(0, 0),
                       I_2 = c(0, 0),
                       S_3  = c(0, 0),
                       I_3  = c(0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         select     = 0,
                         node       = 1,
                         dest       = 0,
                         n          = 1,
                         proportion = 0),
                    .Names = c("event", "time", "select", "node", "dest", "n", "proportion"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(init,
               tspan     = 0:2,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_q1   = 1,
               beta_q2   = 1,
               beta_q3   = 1,
               beta_q4   = 1,
               epsilon   = 0)

result <- run(model, verbose = 0)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L), .Dim = c(12L, 3L))

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(identical(result@U, U))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## Two individuals start in susceptible state in node = 1, with a zero
## probability of becoming infected.
##
## At t = 1, Nkind = 1, one individual is moved to node = 0.
init <- structure(list(id  = c(0, 1),
                       S_1 = c(0, 2),
                       I_1 = c(0, 0),
                       S_2 = c(0, 0),
                       I_2 = c(0, 0),
                       S_3 = c(0, 0),
                       I_3  = c(0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         select     = 0,
                         node       = 1,
                         dest       = 0,
                         n          = 1,
                         proportion = 0),
                    .Names = c("event", "time", "select", "node", "dest", "n", "proportion"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(init,
               tspan     = 0:2,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 0,
               beta_q1   = 1,
               beta_q2   = 1,
               beta_q3   = 1,
               beta_q4   = 1,
               epsilon   = 0)

result <- run(model, verbose = 0)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 1L,
                 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L,
                 0L, 1L, 0L, 0L, 0L, 0L, 0L), .Dim = c(12L, 3L))

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(identical(result@U, U))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))

## 2 Nodes
## 3 Age categories
## 2 Disease-states: Susceptible & Infected
##
## Two individuals start in susceptible state and 8 individuals in
## infected state in node = 1, with a zero probability of becoming
## infected.
##
## At t = 1, Nkind = 2, one individual is moved to node = 0.
init <- structure(list(id  = c(0, 1),
                       S_1 = c(0, 2),
                       I_1 = c(0, 8),
                       S_2 = c(0, 0),
                       I_2 = c(0, 0),
                       S_3 = c(0, 0),
                       I_3 = c(0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -2L), class = "data.frame")

events <- structure(list(event      = 3,
                         time       = 1,
                         select     = 0,
                         node       = 1,
                         dest       = 0,
                         n          = 1,
                         proportion = 0),
                    .Names = c("event", "time", "select", "node", "dest", "n", "proportion"),
                    row.names = c(NA, -1L), class = "data.frame")

model <- SISe3(init,
               tspan     = 0:2,
               events    = events,
               phi       = rep(0, 2),
               upsilon_1 = 0,
               upsilon_2 = 0,
               upsilon_3 = 0,
               gamma_1   = 0,
               gamma_2   = 0,
               gamma_3   = 0,
               alpha     = 0,
               beta_q1   = 1,
               beta_q2   = 1,
               beta_q3   = 1,
               beta_q4   = 1,
               epsilon   = 0)

result <- run(model, verbose = 0, seed = 123L)

U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 2L, 8L, 0L, 0L, 0L, 0L, 0L,
                 1L, 0L, 0L, 0L, 0L, 2L, 7L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L,
                 0L, 2L, 7L, 0L, 0L, 0L, 0L), .Dim = c(12L, 3L))

stopifnot(identical(model@G, result@G))
stopifnot(identical(model@N, result@N))
stopifnot(identical(result@U, U))
stopifnot(identical(model@Nn, result@Nn))
stopifnot(identical(model@data, result@data))
stopifnot(identical(model@sd, result@sd))
stopifnot(identical(model@tspan, result@tspan))
stopifnot(identical(model@u0, result@u0))
stopifnot(identical(model@events, result@events))
