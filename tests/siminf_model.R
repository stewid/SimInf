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
                       data  = matrix(rep(0, Nn), nrow = 1),
                       sd    = rep(0L, Nn),
                       tspan = as.numeric(1),
                       u0    = u0))

tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N,
                       U     = U,
                       Nn    = Nn,
                       data  = matrix(rep(0, Nn), nrow = 1),
                       sd    = rep(0L, Nn),
                       tspan = as.numeric(c(3, 2, 1)),
                       u0    = u0))

## Check u0
tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N,
                       U     = U,
                       Nn    = Nn,
                       data  = matrix(rep(0, Nn), nrow = 1),
                       sd    = rep(0L, Nn),
                       tspan = as.numeric(1:10),
                       u0    = u0 * -1L))

## Check N
tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N * 1.1,
                       U     = U,
                       Nn    = Nn,
                       data  = matrix(rep(0, Nn), nrow = 1),
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
                       data  = matrix(rep(0, Nn), nrow = 1),
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
                       data  = matrix(rep(0, Nn), nrow = 1),
                       sd    = sd,
                       tspan = as.numeric(1:10),
                       u0    = u0))

## Check data
sd <- rep(0L, Nn)

data <- matrix(rep(0, Nn), nrow = 1)
data <- data[, 1:3, drop = FALSE]

## Wrong size of data matrix
tools::assertError(new("siminf_model",
                       G     = G,
                       N     = N,
                       U     = U,
                       Nn    = Nn,
                       data  = data,
                       sd    = sd,
                       tspan = as.numeric(1:10),
                       u0    = u0))

## Check initial state
init <- structure(list(id =      c(0, 1, 2, 3, 4, 5),
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

## Both u0 and init are NULL
tools::assertError(siminf_model())

## Both u0 and init are non NULL
tools::assertError(siminf_model(init = init, u0 = u0))

## Nn must be equal to number of nodes
tools::assertError(siminf_model(u0 = u0, Nn = 7))
