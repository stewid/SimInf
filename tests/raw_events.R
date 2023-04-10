## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2022 Ivana Rodriguez Ewerlöf
## Copyright (C) 2015 -- 2023 Stefan Widgren
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

## Check to pass vectors of different lengths.
res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    integer(0),
    c(0L, 3L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'event' must be an integer vector with length 0.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    c(1L, 1L),
    c(0L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'event' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    c(1L, 1L),
    c(0L, 3L),
    c(1L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'time' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    c(1L, 1L),
    c(0L, 3L),
    c(1L, 2L),
    c(1L),
    c(0L, 2L)))
check_error(res, "'node' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    c(1L, 1L),
    c(0L, 3L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L)))
check_error(res, "'dest' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_raw_events,
    c(1L, 1L),
    c(0L, 2L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'event[2]' is invalid.")

## Check raw events.
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 3L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 1L, 2L, 2L),
    dest  = c(0L, 2L, 2L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(1L, 2L, 4L),
    node  = c(1L, 1L, 2L),
    dest  = c(0L, 2L, 0L))

stopifnot(identical(events_obs, events_exp))

events$id[2] <- NA_integer_
res <- assertError(raw_events(events))
check_error(
    res,
    "'events$id' must be an integer or character vector with non-NA values.")
events$id[2] <- 1L

events$node[2] <- 1.1
res <- assertError(raw_events(events))
check_error(
    res,
    "'events$node' must be an integer or character vector with non-NA values.")
events$node <- c(1L, 1L, 2L, 2L)

events$dest <- as.Date(events$dest, origin = "1970-01-01")
res <- assertError(raw_events(events))
check_error(
    res,
    "'events$dest' must be an integer or character vector with non-NA values.")
events$dest <- c(0L, 2L, 2L, 0L)

## Testing animal with only one enter event, keep
events <- data.frame(
    id    = 1L,
    event = 1L,
    time  = 1L,
    node  = 1L,
    dest  = 0L)

stopifnot(identical(events, as.data.frame(raw_events(events))))

## Testing animal with only one exit event, keep
events <- data.frame(
    id    = 1L,
    event = 0L,
    time  = 1L,
    node  = 1L,
    dest  = 0L)

stopifnot(identical(events, as.data.frame(raw_events(events))))

## Testing animal with only one external transfer event, keep
events <- data.frame(
    id    = 1L,
    event = 3L,
    time  = 1L,
    node  = 1L,
    dest  = 2L)

stopifnot(identical(events, as.data.frame(raw_events(events))))

## Testing animal with two enter events, keep first
events <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 1L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = 1L,
    event = 1L,
    time  = 1L,
    node  = 1L,
    dest  = 0L)

stopifnot(identical(events_obs, events_exp))

## Testing animal with two exit events, keep first
events <- data.frame(
    id    = c(1L, 1L),
    event = c(0L, 0L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = 1L,
    event = 0L,
    time  = 1L,
    node  = 1L,
    dest  = 0L)

stopifnot(identical(events_obs, events_exp))

## Testing animal with two enter events and exit, keep path
events <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 1L, 0L),
    time  = c(1L, 2L, 3L),
    node  = c(1L, 2L, 2L),
    dest  = c(0L, 0L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(2L, 3L),
    node  = c(2L, 2L),
    dest  = c(0L, 0L))

stopifnot(identical(events_obs, events_exp))

## Testing animal with two enter events, a movement and an exit, keep
## path
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 1L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 2L, 1L, 3L),
    dest  = c(0L, 0L, 3L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(1L, 3L, 4L),
    node  = c(1L, 1L, 3L),
    dest  = c(0L, 3L, 0L))

stopifnot(identical(events_obs, events_exp))

## Testing animal with one enter event and two exit events, keep path
events <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 0L, 0L),
    time  = c(1L, 2L, 3L),
    node  = c(1L, 2L, 1L),
    dest  = c(0L, 0L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(1L, 3L),
    node  = c(1L, 1L),
    dest  = c(0L, 0L))

stopifnot(identical(events_obs, events_exp))

## Testing animal with one enter event and two exit events, keep path
events <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 0L, 0L),
    time  = c(1L, 2L, 3L),
    node  = c(1L, 1L, 2L),
    dest  = c(0L, 0L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 0L))

stopifnot(identical(events_obs, events_exp))

## Testing animal with another event after exit event, exit event
## should be last
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 3L, 0L, 3L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 1L, 2L, 2L),
    dest  = c(0L, 2L, 0L, 3L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(1L, 2L, 3L),
    node  = c(1L, 1L, 2L),
    dest  = c(0L, 2L, 0L))

stopifnot(identical(events_obs, events_exp))

## Testing animal with another event after exit event,
## no path to exit, don't keep
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 3L, 0L, 3L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 1L, 3L, 2L),
    dest  = c(0L, 2L, 0L, 3L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = integer(0),
    event = integer(0),
    time  = integer(0),
    node  = integer(0),
    dest  = integer(0))

stopifnot(identical(events_obs, events_exp))

## Testing animal with another event before enter event, enter event
## should be first
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(3L, 1L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 1L, 1L, 2L),
    dest  = c(2L, 0L, 2L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(2L, 3L, 4L),
    node  = c(1L, 1L, 2L),
    dest  = c(0L, 2L, 0L))

stopifnot(identical(events_obs, events_exp))

## Testing animal with another event before enter event, keep path if
## starting on enter event and ending with exit
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(3L, 1L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 2L, 1L, 2L),
    dest  = c(2L, 0L, 2L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(2L, 4L),
    node  = c(2L, 2L),
    dest  = c(0L, 0L))

stopifnot(identical(events_obs, events_exp))

## Testing animal with no path from enter to exit event, don't keep
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(3L, 1L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 2L, 2L, 3L),
    dest  = c(2L, 0L, 1L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = integer(0),
    event = integer(0),
    time  = integer(0),
    node  = integer(0),
    dest  = integer(0))

stopifnot(identical(events_obs, events_exp))

## Testing animal with no enter event, keep path
events <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(3L, 3L, 0L),
    time  = c(1L, 2L, 3L),
    node  = c(1L, 2L, 1L),
    dest  = c(2L, 1L, 0L))

events_obs <- as.data.frame(raw_events(events))

stopifnot(identical(events_obs, events))

## Testing animal with no enter or exit event, keep path
events <- data.frame(
    id    = c(1L, 1L),
    event = c(3L, 3L),
    time  = c(1L, 2L),
    node  = c(1L, 2L),
    dest  = c(2L, 3L))

events_obs <- as.data.frame(raw_events(events))

stopifnot(identical(events_obs, events))

## Testing animal with no exit event, keep path
events <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 3L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 2L))

events_obs <- as.data.frame(raw_events(events))

stopifnot(identical(events_obs, events))

## Testing animal with only enter and exit event, keep path
events <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 0L))

events_obs <- as.data.frame(raw_events(events))

stopifnot(identical(events_obs, events))

## Testing animal with enter and exit event in wrong order, don't keep
events <- data.frame(
    id    = c(1L, 1L),
    event = c(0L, 1L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 0L))

events_obs <- as.data.frame(raw_events(events))

events_exp <- data.frame(
    id    = integer(0),
    event = integer(0),
    time  = integer(0),
    node  = integer(0),
    dest  = integer(0))

stopifnot(identical(events_obs, events_exp))
