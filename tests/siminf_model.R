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
library(Matrix)

## Initialize test data
N <- Matrix(c(-1,  0,  0,
               1,  0,  0,
               0, -1,  0,
               0,  1,  0,
               0,  0, -1,
               0,  0,  1),
            nrow   = 6,
            ncol   = 3,
            byrow  = TRUE,
            sparse = TRUE)

Nn <- 6L

G <- as(Matrix(c(1, 0, 0,
                 0, 1, 0,
                 0, 0, 1),
               nrow   = 3,
               ncol   = 3,
               byrow  = TRUE,
               sparse = TRUE),
        "dgCMatrix")

u0 <- structure(c(0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2,
                  0, 3, 0, 3, 0, 3, 0, 4, 0, 4, 0, 4, 0, 5, 0, 5, 0, 5, 0),
                .Dim = c(6L, 6L))
storage.mode(u0) <- "integer"

U <- matrix(nrow = 0, ncol = 0)
storage.mode(U) <- "integer"

## Check tspan
tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N,
                       U     = U,
                       Nn    = Nn,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       sd    = rep(0L, Nn),
                       tspan = as.numeric(1),
                       u0    = u0))

tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N,
                       U     = U,
                       Nn    = Nn,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       sd    = rep(0L, Nn),
                       tspan = as.numeric(c(3, 2, 1)),
                       u0    = u0))

## Check u0
tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N,
                       U     = U,
                       Nn    = Nn,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       sd    = rep(0L, Nn),
                       tspan = as.numeric(1:10),
                       u0    = u0 * -1L))

## Change storage mode of u0 to double.
## Should not raise error
u0_double <- u0
storage.mode(u0_double) <- "double"
siminf_model(G     = G,
             N     = N,
             U     = U,
             Nn    = Nn,
             ldata = matrix(rep(0, Nn), nrow = 1),
             sd    = rep(0L, Nn),
             tspan = as.numeric(1:10),
             u0    = u0_double)

## Change storage mode of u0 to double and change to non-integer values.
## Should raise error
u0_double <- u0
storage.mode(u0_double) <- "double"
u0_double <- 1.2 * u0_double
tools::assertError(siminf_model(G     = G,
                                N     = N,
                                U     = U,
                                Nn    = Nn,
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                sd    = rep(0L, Nn),
                                tspan = as.numeric(1:10),
                                u0    = u0_double))

## Check N
tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N * 1.1,
                       U     = U,
                       Nn    = Nn,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       sd    = rep(0L, Nn),
                       tspan = as.numeric(1:10),
                       u0    = u0))

## Check G
## Error: Wrong size of dependency graph
tools::assertError(new("siminf_model",
                       G     = G[-1,],
                       N     = N,
                       U     = U,
                       Nn    = Nn,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       sd    = rep(0L, Nn),
                       tspan = as.numeric(1:10),
                       u0    = u0))

## Check sd
sd <- rep(0L, Nn)
sd <- sd[1:3]

## Wrong size of subdomain vector
tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N,
                       U     = U,
                       Nn    = Nn,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       sd    = sd,
                       tspan = as.numeric(1:10),
                       u0    = u0))

## Check gdata
tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N,
                       U     = U,
                       Nn    = Nn,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       gdata = 1L,
                       sd    = rep(0L, Nn),
                       tspan = as.numeric(1:10),
                       u0    = u0))

## Check ldata
sd <- rep(0L, Nn)

ldata <- matrix(rep(0, Nn), nrow = 1)
ldata <- ldata[, 1:3, drop = FALSE]

## Wrong size of ldata matrix
tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N,
                       U     = U,
                       Nn    = Nn,
                       ldata = ldata,
                       sd    = sd,
                       tspan = as.numeric(1:10),
                       u0    = u0))

## Check initial state
init <- structure(list(id  = c(0, 1, 2, 3, 4, 5),
                       S_1 = c(0, 1, 2, 3, 4, 5),
                       I_1 = c(0, 0, 0, 0, 0, 0),
                       S_2 = c(0, 1, 2, 3, 4, 5),
                       I_2 = c(0, 0, 0, 0, 0, 0),
                       S_3 = c(0, 1, 2, 3, 4, 5),
                       I_3 = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -6L), class = "data.frame")

## Both u0 and init are NULL
tools::assertError(siminf_model())

## Both u0 and init are non NULL
tools::assertError(siminf_model(init = init, u0 = u0))

## Nn must be equal to number of nodes
tools::assertError(siminf_model(u0 = u0, Nn = 7))

## Check show method without events
show_expected <- c("Epidemiological model:", "G: 2 x 2", "N: 2 x 2", "U: 0 x 0",
                   "V: 0 x 0", "Nn: 1", "ldata: 0 x 1", "gdata: 1 x 8", "tspan: 1 x 1000",
                   "u0: 2 x 1", "v0: 1 x 1", "",
                   "External events:", "E: 2 x 2", "S: 0 x 0", "event: 0 x 0",
                   "time: 0 x 0", "node: 0 x 0", "dest: 0 x 0", "n: 0 x 0",
                   "proportion: 0 x 0", "select: 0 x 0", "shift: 0 x 0")

show_observed <- capture.output(show(demo_model()))

stopifnot(identical(show_observed, show_expected))

## Check show method with events
init <- structure(list(id  = c(0, 1, 2, 3, 4, 5),
                       S_1 = c(0, 1, 2, 3, 4, 5),
                       I_1 = c(0, 0, 0, 0, 0, 0),
                       S_2 = c(0, 1, 2, 3, 4, 5),
                       I_2 = c(0, 0, 0, 0, 0, 0),
                       S_3 = c(0, 1, 2, 3, 4, 5),
                       I_3 = c(0, 0, 0, 0, 0, 0)),
                  .Names = c("id", "S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                  row.names = c(NA, -6L),
                  class = "data.frame")

events <- structure(list(event      = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         node       = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest       = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select     = c(3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5),
                         shift      = c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(init,
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
               beta_q1   = 1,
               beta_q2   = 1,
               beta_q3   = 1,
               beta_q4   = 1,
               epsilon   = 1)

show_expected <- c("Epidemiological model:", "G: 6 x 6", "N: 6 x 6", "U: 0 x 0",
                   "V: 0 x 0", "Nn: 6", "ldata: 0 x 6", "gdata: 1 x 12",
                   "tspan: 1 x 11", "u0: 6 x 6", "v0: 1 x 6", "",
                   "External events:", "E: 6 x 6", "S: 6 x 2", "event: 1 x 15",
                   "time: 1 x 15", "node: 1 x 15", "dest: 1 x 15", "n: 1 x 15",
                   "proportion: 1 x 15", "select: 1 x 15", "shift: 1 x 15")

show_observed <- capture.output(show(model))

stopifnot(identical(show_observed, show_expected))

## Check U. Change storage mode of U to double.
## Should not raise error
U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L, 2L, 0L, 1L, 1L,
                 1L, 1L, 2L, 1L, 3L, 0L, 2L, 1L, 2L, 2L, 0L, 4L, 1L, 3L, 2L, 3L,
                 3L, 2L, 1L, 4L, 6L, 9L, 7L, 8L, 4L, 11L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 7L, 8L, 7L, 8L, 5L, 10L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 8L, 7L, 6L, 9L,
                 8L, 7L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 5L, 10L, 4L, 11L, 6L, 9L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 7L, 8L, 5L, 10L, 7L, 8L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 11L, 5L, 10L, 3L, 12L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 8L, 7L, 5L,
                 10L, 4L, 11L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 6L, 9L, 2L, 13L, 4L, 11L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 5L, 10L, 2L, 13L, 7L, 8L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 9L, 6L, 2L, 13L, 6L, 9L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
               .Dim = c(36L, 11L))

U_double <- U
storage.mode(U_double) <- "double"

siminf_model(G     = G,
             N     = N,
             U     = U_double,
             Nn    = Nn,
             ldata = matrix(rep(0, Nn), nrow = 1),
             sd    = rep(0L, Nn),
             tspan = as.numeric(1:10),
             u0    = u0)

## Check U. Change storage mode of U to double and change to non-integer values.
## Should raise error
U_double <- U
storage.mode(U_double) <- "double"
U_double <- U_double * 1.2

tools::assertError(siminf_model(G     = G,
                                N     = N,
                                U     = U_double,
                                Nn    = Nn,
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                sd    = rep(0L, Nn),
                                tspan = as.numeric(1:10),
                                u0    = u0))

## Check U. Should not raise an error if U is an integer vector of length 0
siminf_model(G     = G,
             N     = N,
             U     = integer(0),
             Nn    = Nn,
             ldata = matrix(rep(0, Nn), nrow = 1),
             sd    = rep(0L, Nn),
             tspan = as.numeric(1:10),
             u0    = u0)

## Check U. Should raise error if U is an integer vector of length > 0
tools::assertError(siminf_model(G     = G,
                                N     = N,
                                U     = c(1L),
                                Nn    = Nn,
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                sd    = rep(0L, Nn),
                                tspan = as.numeric(1:10),
                                u0    = u0))

## Check V. Change storage mode of V to double.
## Should not raise error
V <- structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1),
               .Dim = c(6L, 10L))

V_integer <- V
storage.mode(V) <- "integer"

siminf_model(G     = G,
             N     = N,
             U     = U,
             V     = V_integer,
             Nn    = Nn,
             ldata = matrix(rep(0, Nn), nrow = 1),
             sd    = rep(0L, Nn),
             tspan = as.numeric(1:10),
             u0    = u0)

## Check V. Change storage mode of V to character
## Should raise error
V_character <- V
storage.mode(V_character) <- "character"

tools::assertError(siminf_model(G     = G,
                                N     = N,
                                U     = U,
                                V     = V_character,
                                Nn    = Nn,
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                sd    = rep(0L, Nn),
                                tspan = as.numeric(1:10),
                                u0    = u0))

## Check V. Should raise error if V is a vector of length > 0
tools::assertError(siminf_model(G     = G,
                                N     = N,
                                U     = U,
                                V     = 1,
                                Nn    = Nn,
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                sd    = rep(0L, Nn),
                                tspan = as.numeric(1:10),
                                u0    = u0))

## Check V. Should not raise an error if V is an integer vector of length 0
siminf_model(G     = G,
             N     = N,
             U     = U,
             V     = integer(0),
             Nn    = Nn,
             ldata = matrix(rep(0, Nn), nrow = 1),
             sd    = rep(0L, Nn),
             tspan = as.numeric(1:10),
             u0    = u0)
