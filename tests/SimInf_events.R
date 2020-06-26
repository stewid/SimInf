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
library(Matrix)
library(tools)
source("util/check.R")

## For debugging
sessionInfo()

E <- Matrix(c(1, 0, 0, 1, 0, 0,
              0, 0, 0, 1, 0, 0,
              0, 1, 0, 0, 1, 0,
              0, 0, 0, 0, 1, 0,
              0, 0, 1, 0, 0, 1,
              0, 0, 0, 0, 0, 1),
            nrow   = 6,
            ncol   = 6,
            byrow  = TRUE,
            sparse = TRUE,
            dimnames = list(c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                            c("1", "2", "3", "4", "5", "6")))

N <- matrix(c(2, 0,
              2, 0,
              0, 2,
              0, 2,
              0, 0,
              0, 0),
            nrow   = 6,
            ncol   = 2,
            byrow  = TRUE,
            dimnames = list(c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                            c("1", "2")))

## Check valid_SimInf_events_object
events <- SimInf_events(E = E, N = N)
stopifnot(isTRUE(SimInf:::valid_SimInf_events_object(events)))

rownames(events@E) <- tolower(rownames(events@E))
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "'E' and 'N' must have identical compartments."))
rownames(events@E) <- NULL
stopifnot(identical(
    SimInf:::valid_SimInf_events_object(events),
    "'E' and 'N' must have rownames matching the compartments."))

events <- SimInf_events(E = E, N = N)
events@E[1, 1] <- -1
stopifnot(identical(
    SimInf:::valid_SimInf_events_object(events),
    "Select matrix 'E' has negative elements."))

events <- SimInf_events(E = E, N = N)
events@event <- 3L
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "All scheduled events must have equal length."))

events@event <- 3L
events@time <- 1L
events@node <- 2L
events@dest <- 1L
events@n <- 1L
events@proportion <- 0
events@select <- 1L
events@shift <- 1L

events@time <- 0L
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "time must be greater than 0."))
events@time <- 1L

events@event <- -1L
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "event must be in the range 0 <= event <= 3."))
events@event <- 4L
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "event must be in the range 0 <= event <= 3."))
events@event <- 3L

events@node <- 0L
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "'node' must be greater or equal to 1."))
events@node <- 2L

events@dest <- 0L
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "'dest' must be greater or equal to 1."))
events@dest <- 1L

events@proportion <- -1
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "prop must be in the range 0 <= prop <= 1."))
events@proportion <- 2
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "prop must be in the range 0 <= prop <= 1."))
events@proportion <- 0

events@select <- 0L
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "select must be in the range 1 <= select <= Nselect."))
events@select <- 7L
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "select must be in the range 1 <= select <= Nselect."))
events@select <- 1L

events@event <- 2L
events@shift <- 0L
stopifnot(identical(SimInf:::valid_SimInf_events_object(events),
                    "'shift' must be greater or equal to 1."))

## Check that an error is raised when E is NULL and events is
## non-NULL.
events <- data.frame(event = 3, time = 1, node = 2, dest = 1, n = 1,
                     proportion = 1, select = 0, shift = 1)
res <- assertError(SimInf_events(events = events))
check_error(res, "events is not NULL when E is NULL.")

## Check that an error is raised when N is not an integer matrix.
res <- assertError(SimInf_events(N = c("a", "b")))
check_error(res, "'N' must be an integer matrix.")

## Check that an error is raised when N is not an integer matrix.
res <- assertError(SimInf_events(N = matrix(c(1.5, 1.5), nrow = 2)))
check_error(res, "'N' must be an integer matrix.")

## Check that an error is raised when events is not a data.frame.
res <- assertError(SimInf_events(E = E, events = c("a", "b")))
check_error(res, "events must be a data.frame.")

## Check that an error is raised when events contains a non-numeric
## value.
events <- data.frame(event = 3, time = 1, node = 2, dest = 1, n = 1,
                     proportion = "1", select = 0, shift = 1)
res <- assertError(SimInf_events(E = E, events = events))
check_error(res, "Columns in events must be numeric.")

## Check missing columns in events
## Iterate over each column and rename it
lapply(seq_len(8), function(i) {
    events <- data.frame(
        event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
        time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
        dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
        proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
        shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))

    colnames(events)[i] <- "test"
    res <- assertError(SimInf_events(E = E, N = N, events = events))
    check_error(res, "Missing columns in events.")
})

## Check events$event not equal to whole number
events <- data.frame(
    event = c(3.1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
str(events)
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(res, "Columns in events must be integer.")

## Check events$time not equal to whole number
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(1.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
str(events)
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(res, "Columns in events must be integer.")

## Check missing t0 when events$time is a Date vector
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = structure(c(17168, 17169, 17170, 17171, 17172,
                       17173, 17174, 17175, 17176, 17177,
                       17178, 17179, 17180, 17181, 17182),
                     class = "Date"),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(res, "Missing 't0'.")

## Check invalid t0 (length != 1) when events$time is a Date vector
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = structure(c(17168, 17169, 17170, 17171, 17172,
                       17173, 17174, 17175, 17176, 17177,
                       17178, 17179, 17180, 17181, 17182),
                     class = "Date"),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
res <- assertError(SimInf_events(E = E, N = N, events = events, t0 = c(1, 1)))
check_error(res, "Invalid 't0'.")

## Check invalid t0 (!numeric) when events$time is a Date vector
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = structure(c(17168, 17169, 17170, 17171, 17172,
                       17173, 17174, 17175, 17176, 17177,
                       17178, 17179, 17180, 17181, 17182),
                     class = "Date"),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
res <- assertError(SimInf_events(E = E, N = N, events = events, t0 = "1"))
check_error(res, "Invalid 't0'.")

## Check invalid t0 (!NULL) when events$time is not a Date vector
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
res <- assertError(SimInf_events(E = E, N = N, events = events, t0 = 0))
check_error(res, "Invalid 't0'.")

## Check events$time equal to a Date vector
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = structure(c(17168, 17169, 17170, 17171, 17172,
                       17173, 17174, 17175, 17176, 17177,
                       17178, 17179, 17180, 17181, 17182),
                     class = "Date"),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
res <- SimInf_events(E = E, N = N, events = events, t0 = 17166)
stopifnot(identical(res@time, 2:16))

## Check events$time equal to an integer vector
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(17168, 17169, 17170, 17171, 17172,
             17173, 17174, 17175, 17176, 17177,
             17178, 17179, 17180, 17181, 17182),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
res <- SimInf_events(E = E, N = N, events = events)
stopifnot(identical(is.null(names(res@time)), TRUE))

## Check events$select not equal to whole number
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(0.1, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
str(events)
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(res, "Columns in events must be integer.")

## Check events$node not equal to whole number
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(2.2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
str(events)
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(res, "Columns in events must be integer.")

## Check events$node less than one
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
str(events)
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(res, "'node' must be greater or equal to 1", FALSE)

## Check events$dest not equal to whole number
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
str(events)
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(res, "Columns in events must be integer.")

## Check events$dest less than 1
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
str(events)
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(res, "'dest' must be greater or equal to 1", FALSE)

## Check events$n not equal to whole number
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1.1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
str(events)
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(res, "Columns in events must be integer.")

## Check events$event equal to character
events <- data.frame(
    event = c("3", "3", "3", "3", "3", "3", "3", "3", "3",
              "3", "3", "3", "3", "3", "3"),
    time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    shift = c(1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
str(events)
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(
    res,
    "'event' type must be 'enter', 'exit', 'extTrans' or 'intTrans'.")

## Check events$shift not equal to whole number
events <- data.frame(
    event = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2),
    shift = c(1.1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
str(events)
res <- assertError(SimInf_events(E = E, N = N, events = events))
check_error(res, "Columns in events must be integer.")

## Check E and events equal to NULL (default).
events <- new("SimInf_events",
              E = new("dgCMatrix",
                      i = integer(0),
                      p = 0L,
                      Dim = c(0L, 0L),
                      Dimnames = list(NULL, NULL),
                      x = numeric(0),
                      factors = list()),
              N = matrix(integer(0), nrow = 0, ncol = 0),
              event = integer(0),
              time = integer(0),
              node = integer(0),
              dest = integer(0),
              n = integer(0),
              proportion = numeric(0),
              select = integer(0),
              shift = integer(0))
str(events)
stopifnot(identical(SimInf_events(), events))

## Check the SimInf_events plot method. Reduce the run-time by only
## using one year of data
data(events_SISe3)
events_SISe3 <- events_SISe3[events_SISe3$time < 366, ]
events <- SimInf_events(E = E, N = N, events = events_SISe3)
stopifnot(identical(events, show(events)))
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(events)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Test summary method with scheduled events
summary_expected <-
    c("Number of scheduled events: 195919",
      " - Exit: 45562 (n: min = 1 max = 1 avg = 1.0)",
      " - Enter: 45603 (n: min = 1 max = 1 avg = 1.0)",
      " - Internal transfer: 79225 (n: min = 1 max = 4 avg = 1.1)",
      " - External transfer: 25529 (n: min = 1 max = 1 avg = 1.0)")
summary_observed <- capture.output(summary(events))
stopifnot(identical(summary_observed, summary_expected))

## Test summary method with no scheduled events
events <- SimInf_events(E = E, N = N)
summary_expected <-
    c("Number of scheduled events: 0",
      " - Exit: 0",
      " - Enter: 0",
      " - Internal transfer: 0",
      " - External transfer: 0")
summary_observed <- capture.output(summary(events))
stopifnot(identical(summary_observed, summary_expected))

## Check if converting the events to a data.frame results in the same
## as the events data submitted to the SimInf_events function.
events <- data.frame(
    event = c(3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
              3L, 3L, 3L, 3L, 3L, 3L),
    time = c(17168L, 17169L, 17170L, 17171L, 17172L,
             17173L, 17174L, 17175L, 17176L, 17177L,
             17178L, 17179L, 17180L, 17181L, 17182L),
    node = c(2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L,
             5L, 5L, 5L, 6L, 6L, 6L),
    dest = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
             1L, 1L, 1L, 1L, 1L, 1L),
    n = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L,
          4L, 4L, 5L, 5L, 5L),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1),
    select = c(1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L,
               2L, 1L, 1L, 2L, 1L, 1L, 2L),
    shift = c(1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L,
              1L, 2L, 0L, 1L, 2L, 0L))
res <- SimInf_events(E = E, N = N, events = events)
stopifnot(identical(as(res, "data.frame"), events))
stopifnot(identical(as.data.frame(res), events))

## Check that it fails when dest is out of bounds.
u0 <- data.frame(S = c(10, 10), I = c(0, 0), R = c(0, 0))
events <- data.frame(event = 3, time = 2, node = 1, dest = 3,
                     n = 1, proportion = 0, select = 2, shift = 0)
model <- SIR(u0, tspan = seq_len(3), events = events, beta = 0.16,
             gamma = 0.077)
res <- assertError(run(model))
check_error(res, "'dest' is out of bounds.")

## Check that it fails when node is out of bounds.
u0 <- data.frame(S = c(10, 10), I = c(0, 0), R = c(0, 0))
events <- data.frame(event = 3, time = 2, node = 1, dest = 2,
                     n = 1, proportion = 0, select = 2, shift = 0)
model <- SIR(u0, tspan = seq_len(3), events = events, beta = 0.16,
             gamma = 0.077)
model@events@node <- -1L
res <- assertError(.Call(SimInf:::SIR_run, model, NULL))
check_error(res, "'node' is out of bounds.")

## Check that it fails for an invalid event type.
u0 <- data.frame(S = c(10, 10), I = c(0, 0), R = c(0, 0))
events <- data.frame(event = 0, time = 2, node = 1, dest = 0,
                     n = 1, proportion = 0, select = 2, shift = 0)
model <- SIR(u0, tspan = seq_len(3), events = events, beta = 0.16,
             gamma = 0.077)
model@events@event <- 4L
res <- assertError(.Call(SimInf:::SIR_run, model, NULL))
check_error(res, "Undefined event type.")

## Check get/set select_matrix
model <- SIR(cbind(S = 100, I = 10, R = 0), tspan = 1:10, beta = 1,
             gamma = 1)

## Set the select matrix
select_matrix(model) <- matrix(c(1, 0, 0, 1, 1, 1, 0, 0, 1), nrow = 3)

E_expected <- new("dgCMatrix",
                  i = c(0L, 0L, 1L, 2L, 2L),
                  p = c(0L, 1L, 4L, 5L),
                  Dim = c(3L, 3L),
                  Dimnames = list(c("S", "I", "R"), c("1", "2", "3")),
                  x = c(1, 1, 1, 1, 1),
                  factors = list())

## Extract the select matrix from the model
E_observed <- select_matrix(model)

stopifnot(identical(E_expected, E_observed))

m <- matrix(c(1, 0, 0, 1, 1, 1, 0, 0), nrow = 2)
res <- assertError(select_matrix(model) <- m)
check_error(res, "'value' must have one row for each compartment in the model.")

## Check get/set shift_matrix
model <- SIR(cbind(S = 100, I = 10, R = 0), tspan = 1:10, beta = 1, gamma = 1)

## Set the shift matrix
shift_matrix(model) <- matrix(c(2, 1, 0), nrow = 3)

N_expected <- structure(c(2L, 1L, 0L), .Dim = c(3L, 1L),
                        .Dimnames = list(c("S", "I", "R"), "1"))

## Extract the shift matrix from the model
N_observed <- shift_matrix(model)

stopifnot(identical(N_expected, N_observed))

shift_matrix(model) <- NULL
stopifnot(identical(shift_matrix(model),
                    matrix(integer(0), nrow = 0, ncol = 0)))

m <- matrix(c(1, 0), nrow = 2, dimnames = list(c("S", "I")))
res <- assertError(select_matrix(model) <- m)
check_error(res, "'value' must have one row for each compartment in the model.")

m <- matrix(c("1", "0", "0"), nrow = 3)
res <- assertError(shift_matrix(model) <- m)
check_error(res, "'N' must be an integer matrix.")

m <- matrix(c(1.3, 0, 0), nrow = 3)
res <- assertError(shift_matrix(model) <- m)
check_error(res, "'N' must be an integer matrix.")
