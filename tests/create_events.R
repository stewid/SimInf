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
    .Call(siminf:::siminf_external_events, nodes, 10L, rep(0.1, 10), rep(0.1, 10), NULL))
stopifnot(length(grep("Invalid 'nodes' argument",
                      res[[1]]$message)) > 0)

nodes <- integer(0)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, nodes, 10L, rep(0.1, 10), rep(0.1, 10), NULL))
stopifnot(length(grep("Invalid 'nodes' argument",
                      res[[1]]$message)) > 0)

nodes <- 5
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, nodes, 10L, rep(0.1, 10), rep(0.1, 10), NULL))
stopifnot(length(grep("Invalid 'nodes' argument",
                      res[[1]]$message)) > 0)

nodes <- NA_integer_
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, nodes, 10L, rep(0.1, 10), rep(0.1, 10), NULL))
stopifnot(length(grep("Invalid 'nodes' argument",
                      res[[1]]$message)) > 0)

nodes <- c(5L, 10L)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, nodes, 10L, rep(0.1, 10), rep(0.1, 10), NULL))
stopifnot(length(grep("Invalid 'nodes' argument",
                      res[[1]]$message)) > 0)

## Check 'days' argument to C function 'siminf_external_events'
days <- NULL
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, days, 0.1, 0.1, NULL))
stopifnot(length(grep("Invalid 'days' argument",
                      res[[1]]$message)) > 0)

days <- integer()
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, days, 0.1, 0.1, NULL))
stopifnot(length(grep("Invalid 'days' argument",
                      res[[1]]$message)) > 0)

days <- 5
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, days, 0.1, 0.1, NULL))
stopifnot(length(grep("Invalid 'days' argument",
                      res[[1]]$message)) > 0)

days <- NA_integer_
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, days, 0.1, 0.1, NULL))
stopifnot(length(grep("Invalid 'days' argument",
                      res[[1]]$message)) > 0)

days <- c(5L, 10L)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, days, 0.1, 0.1, NULL))
stopifnot(length(grep("Invalid 'days' argument",
                      res[[1]]$message)) > 0)

## Check 'p_edge' argument to C function 'siminf_external_events'
p_edge <- NULL
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5), NULL))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

p_edge <- numeric(0)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5), NULL))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

p_edge <- 5
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5), NULL))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

p_edge <- NA_real_
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5), NULL))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

p_edge <- c(rep(0.1, 4), NA_real_)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5), NULL))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

p_edge <- c(rep(0.1, 4), Inf)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5), NULL))
stopifnot(length(grep("Invalid 'p_edge' argument",
                      res[[1]]$message)) > 0)

p_edge <- rep(0.0, 5)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5), NULL))
stopifnot(length(grep("Invalid 'p_edge': Must be in interval 0 < p_edge < 1",
                      res[[1]]$message)) > 0)

p_edge <- rep(1.0, 5)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, p_edge, rep(0.1, 5), NULL))
stopifnot(length(grep("Invalid 'p_edge': Must be in interval 0 < p_edge < 1",
                      res[[1]]$message)) > 0)

## Check 'mu' argument to C function 'siminf_external_events'
mu <- NULL
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu, NULL))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)

mu <- numeric(0)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu, NULL))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)

mu <- 5
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu, NULL))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)

mu <- NA_real_
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu, NULL))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)

mu <- c(rep(0.1, 4), NA_real_)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu, NULL))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)

mu <- c(rep(0.1, 4), Inf)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 5L, rep(0.1, 5), mu, NULL))
stopifnot(length(grep("Invalid 'mu' argument",
                      res[[1]]$message)) > 0)

## Check 'seed' argument to C function 'siminf_external_events'
seed <- integer(0)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 10L, rep(0.1, 10), rep(0.1, 10), seed))
stopifnot(length(grep("Invalid 'seed' argument",
                      res[[1]]$message)) > 0)

seed <- 5
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 10L, rep(0.1, 10), rep(0.1, 10), seed))
stopifnot(length(grep("Invalid 'seed' argument",
                      res[[1]]$message)) > 0)

seed <- NA_integer_
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 10L, rep(0.1, 10), rep(0.1, 10), seed))
stopifnot(length(grep("Invalid 'seed' argument",
                      res[[1]]$message)) > 0)

seed <- c(5L, 10L)
res <- tools::assertError(
    .Call(siminf:::siminf_external_events, 10L, 10L, rep(0.1, 10), rep(0.1, 10), seed))
stopifnot(length(grep("Invalid 'seed' argument",
                      res[[1]]$message)) > 0)

## Create events
events_exp <- structure(list(
    event = c(3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
              3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
              3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L),
    time = c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L,
             3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 5L, 5L, 5L, 6L, 6L, 6L,
             6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L, 8L, 9L, 9L, 9L, 9L, 9L),
    node = c(3L, 6L, 7L, 8L, 8L, 0L, 0L, 6L, 7L, 8L, 0L, 0L, 2L, 8L, 8L, 0L,
             4L, 5L, 5L, 6L, 7L, 8L, 9L, 0L, 2L, 6L, 0L, 5L, 7L, 0L, 0L, 6L,
             7L, 7L, 8L, 0L, 4L, 6L, 8L, 9L, 0L, 0L, 0L, 0L, 7L, 8L),
    dest = c(8L, 8L, 6L, 5L, 9L, 1L, 9L, 7L, 9L, 9L, 5L, 6L, 9L, 3L, 9L, 7L,
             7L, 7L, 8L, 7L, 5L, 4L, 8L, 3L, 9L, 7L, 7L, 9L, 6L, 3L, 8L, 7L,
             4L, 9L, 6L, 6L, 8L, 7L, 3L, 4L, 8L, 5L, 6L, 9L, 9L, 9L),
    n = c(8L, 4L, 5L, 3L, 4L, 7L, 4L, 4L, 4L, 3L, 8L, 4L, 4L, 6L, 5L, 6L, 7L,
          4L, 3L, 2L, 6L, 5L, 4L, 5L, 4L, 3L, 1L, 3L, 4L, 2L, 2L, 4L, 6L, 6L,
          3L, 1L, 7L, 6L, 6L, 4L, 5L, 2L, 1L, 4L, 3L, 3L),
    select = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
               0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
               0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
    shift = c(-1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L,
              -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L,
              -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L, -1L,
              -1L, -1L, -1L, -1L, -1L, -1L, -1L)),
    .Names = c("event", "time", "node", "dest", "n", "select", "shift"))

events_obs <- .Call(siminf:::siminf_external_events, 10L, 10L, rep(0.1, 10), rep(4.25, 10), 123L)
stopifnot(identical(events_obs, events_exp))
