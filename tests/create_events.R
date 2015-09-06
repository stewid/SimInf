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

## Check 'nodes' argument to C function 'siminf_external_events'
nodes <- NULL
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, nodes, 10, 0.1, 0.1))
stopifnot(length(grep("Invalid 'nodes' argument",
                      res[[1]]$message)) > 0)

nodes <- integer(0)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, nodes, 10, 0.1, 0.1))
stopifnot(length(grep("Invalid 'nodes' argument",
                      res[[1]]$message)) > 0)

nodes <- 5
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, nodes, 10, 0.1, 0.1))
stopifnot(length(grep("Invalid 'nodes' argument",
                      res[[1]]$message)) > 0)

nodes <- NA_integer_
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, nodes, 10, 0.1, 0.1))
stopifnot(length(grep("Invalid 'nodes' argument",
                      res[[1]]$message)) > 0)

nodes <- c(5L, 10L)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, nodes, 10, 0.1, 0.1))
stopifnot(length(grep("Invalid 'nodes' argument",
                      res[[1]]$message)) > 0)

## Check 'days' argument to C function 'siminf_external_events'
days <- NULL
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, days, 0.1, 0.1))
stopifnot(length(grep("Invalid 'days' argument",
                      res[[1]]$message)) > 0)

days <- integer()
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, days, 0.1, 0.1))
stopifnot(length(grep("Invalid 'days' argument",
                      res[[1]]$message)) > 0)

days <- 5
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, days, 0.1, 0.1))
stopifnot(length(grep("Invalid 'days' argument",
                      res[[1]]$message)) > 0)

days <- NA_integer_
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, days, 0.1, 0.1))
stopifnot(length(grep("Invalid 'days' argument",
                      res[[1]]$message)) > 0)

days <- c(5L, 10L)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, days, 0.1, 0.1))
stopifnot(length(grep("Invalid 'days' argument",
                      res[[1]]$message)) > 0)

## Check 'p_edge' argument to C function 'siminf_external_events'
p_edge <- NULL
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5)))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

p_edge <- numeric(0)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5)))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

p_edge <- 5
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5)))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

p_edge <- NA_real_
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5)))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

p_edge <- c(0.1, NA_real_)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5)))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

## Check 'mu' argument to C function 'siminf_external_events'
mu <- NULL
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)

mu <- numeric(0)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)

mu <- 5
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)

mu <- NA_real_
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)

mu <- c(0.1, NA_real_)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)
