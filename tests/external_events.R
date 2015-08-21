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

E <- Matrix(c(1, 0, 0, 1, 0, 0,
              0, 0, 0, 1, 0, 0,
              0, 1, 0, 0, 1, 0,
              0, 0, 0, 0, 1, 0,
              0, 0, 1, 0, 0, 1,
              0, 0, 0, 0, 0, 1),
            nrow   = 6,
            ncol   = 6,
            byrow  = TRUE,
            sparse = TRUE)

S <- Matrix(c(2, 0,
              2, 0,
              0, 2,
              0, 2,
              0, 0,
              0, 0),
            nrow   = 6,
            ncol   = 2,
            byrow  = TRUE,
            sparse = TRUE)

## Check events$event not equal to whole number
events <- structure(list(event = c(3.1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         node = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         shift = c(0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1)),
                    .Names = c("event", "time", "node", "dest", "n", "proportion", "select", "shift"),
                    row.names = c(NA, -15L), class = "data.frame")
str(events)
res <- tools::assertError(external_events(E      = E,
                                          S      = S,
                                          events = events))
stopifnot(length(grep("Columns in events must be integer",
                      res[[1]]$message)) > 0)

## Check events$time not equal to whole number
events <- structure(list(event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         time = c(1.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         node = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         shift = c(0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1)),
                    .Names = c("event", "time", "node", "dest", "n", "proportion", "select", "shift"),
                    row.names = c(NA, -15L), class = "data.frame")
str(events)
res <- tools::assertError(external_events(E      = E,
                                          S      = S,
                                          events = events))
stopifnot(length(grep("Columns in events must be integer",
                      res[[1]]$message)) > 0)

## Check events$select not equal to whole number
events <- structure(list(event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         node = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0.1, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         shift = c(0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1)),
                    .Names = c("event", "time", "node", "dest", "n", "proportion", "select", "shift"),
                    row.names = c(NA, -15L), class = "data.frame")
str(events)
res <- tools::assertError(external_events(E      = E,
                                          S      = S,
                                          events = events))
stopifnot(length(grep("Columns in events must be integer",
                      res[[1]]$message)) > 0)

## Check events$node not equal to whole number
events <- structure(list(event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         node = c(1.1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         shift = c(0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1)),
                    .Names = c("event", "time", "node", "dest", "n", "proportion", "select", "shift"),
                    row.names = c(NA, -15L), class = "data.frame")
str(events)
res <- tools::assertError(external_events(E      = E,
                                          S      = S,
                                          events = events))
stopifnot(length(grep("Columns in events must be integer",
                      res[[1]]$message)) > 0)

## Check events$dest not equal to whole number
events <- structure(list(event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         node = c(1.1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest = c(0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         shift = c(0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1)),
                    .Names = c("event", "time", "node", "dest", "n", "proportion", "select", "shift"),
                    row.names = c(NA, -15L), class = "data.frame")
str(events)
res <- tools::assertError(external_events(E      = E,
                                          S      = S,
                                          events = events))
stopifnot(length(grep("Columns in events must be integer",
                      res[[1]]$message)) > 0)

## Check events$n not equal to whole number
events <- structure(list(event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         node = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n = c(1.1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         shift = c(0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1)),
                    .Names = c("event", "time", "node", "dest", "n", "proportion", "select", "shift"),
                    row.names = c(NA, -15L), class = "data.frame")
str(events)
res <- tools::assertError(external_events(E      = E,
                                          S      = S,
                                          events = events))
stopifnot(length(grep("Columns in events must be integer",
                      res[[1]]$message)) > 0)

## Check events$event equal to character
events <- structure(list(event = c("3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3"),
                         time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         node = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         shift = c(0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1)),
                    .Names = c("event", "time", "node", "dest", "n", "proportion", "select", "shift"),
                    row.names = c(NA, -15L), class = "data.frame")
str(events)
res <- tools::assertError(external_events(E      = E,
                                          S      = S,
                                          events = events))
stopifnot(length(grep("Columns in events must be numeric",
                      res[[1]]$message)) > 0)

## Check events$shift not equal to whole number
events <- structure(list(event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                         time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         node = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         dest = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                         n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
                         proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
                         shift = c(0.1, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1)),
                    .Names = c("event", "time", "node", "dest", "n", "proportion", "select", "shift"),
                    row.names = c(NA, -15L), class = "data.frame")
str(events)
res <- tools::assertError(external_events(E      = E,
                                          S      = S,
                                          events = events))
stopifnot(length(grep("Columns in events must be integer",
                      res[[1]]$message)) > 0)

## Check E and events equal to NULL (default).
events <- new("external_events",
              E = new("dgCMatrix",
                  i = integer(0),
                  p = 0L,
                  Dim = c(0L, 0L),
                  Dimnames = list(NULL, NULL),
                  x = numeric(0),
                  factors = list()),
              S = new("dgCMatrix",
                  i = integer(0),
                  p = 0L,
                  Dim = c(0L, 0L),
                  Dimnames = list(NULL, NULL),
                  x = numeric(0),
                  factors = list()),
              event = integer(0),
              time = integer(0),
              node = integer(0),
              dest = integer(0),
              n = integer(0),
              proportion = numeric(0),
              select = integer(0),
              shift = integer(0))
str(events)
stopifnot(identical(external_events(), events))
