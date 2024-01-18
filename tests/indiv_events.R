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

## Specify the number of threads to use.
set_num_threads(1)

## For debugging
sessionInfo()

## Check to pass vectors of different lengths.
res <- assertError(.Call(
    SimInf:::SimInf_clean_indiv_events,
    integer(0),
    c(0L, 3L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'event' must be an integer vector with length 0.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_indiv_events,
    c(1L, 1L),
    c(0L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'event' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_indiv_events,
    c(1L, 1L),
    c(0L, 3L),
    c(1L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'time' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_indiv_events,
    c(1L, 1L),
    c(0L, 3L),
    c(1L, 2L),
    c(1L),
    c(0L, 2L)))
check_error(res, "'node' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_indiv_events,
    c(1L, 1L),
    c(0L, 3L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L)))
check_error(res, "'dest' must be an integer vector with length 2.")

res <- assertError(.Call(
    SimInf:::SimInf_clean_indiv_events,
    c(1L, 1L),
    c(0L, 2L),
    c(1L, 2L),
    c(1L, 1L),
    c(0L, 2L)))
check_error(res, "'event[2]' is invalid.")

## Check individual events.
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 3L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 1L, 2L, 2L),
    dest  = c(0L, 2L, 2L, 1L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(1L, 2L, 4L),
    node  = c(1L, 1L, 2L),
    dest  = c(NA_integer_, 2L, NA_integer_))

stopifnot(identical(events_obs, events_exp))

events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 3L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c("1", "1", "2", "2"),
    dest  = c("0", "2", "2", "1"))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(1L, 2L, 4L),
    node  = c("1", "1", "2"),
    dest  = c(NA_character_, "2", NA_character_))

stopifnot(identical(events_obs, events_exp))

events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 3L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 1L, 2L, 2L),
    dest  = c(0L, 2L, 2L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(1L, 2L, 4L),
    node  = c(1L, 1L, 2L),
    dest  = c(NA_integer_, 2L, NA_integer_))

stopifnot(identical(events_obs, events_exp))

events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 3L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c("A", "A", "B", "B"),
    dest  = c("0", "B", "B", "0"))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(1L, 2L, 4L),
    node  = c("A", "A", "B"),
    dest  = c(NA_character_, "B", NA_character_))

stopifnot(identical(events_obs, events_exp))

events <- data.frame(
    id    = c("A", "A", "A", "A"),
    event = c("enter", "extTrans", "extTrans", "exit"),
    time  = c("2019-02-02", "2020-03-07", "2021-04-14", "2022-05-11"),
    node  = c(1L, 1L, 2L, 2L),
    dest  = c(0L, 2L, 2L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id = c("A", "A", "A"),
    event = c("enter", "extTrans", "exit"),
    time  = as.Date(c("2019-02-02", "2020-03-07", "2022-05-11")),
    node = c(1L, 1L, 2L),
    dest = c(NA_integer_, 2L, NA_integer_))

stopifnot(identical(events_obs, events_exp))

events$id[2] <- NA_integer_
res <- assertError(individual_events(events))
check_error(
    res,
    "'id' must be an integer or character vector with non-NA values.")
events$id[2] <- 1L

events$node[2] <- 1.1
res <- assertError(individual_events(events))
check_error(
    res,
    "'node' and 'dest' must both be integer or character.")
events$node <- c(1L, 1L, 2L, 2L)

events$dest <- as.Date(events$dest, origin = "1970-01-01")
res <- assertError(individual_events(events))
check_error(
    res,
    "'node' and 'dest' must both be integer or character.")
events$dest <- c(0L, 2L, 2L, 0L)

events <- data.frame(
    id    = c("A", "A", "A", "A"),
    event = c("enter", "extTrans", "extTrans", "exit"),
    time  = c("2001-02-01", 2L, 3L, 4L),
    node  = c(1L, 1L, 2L, 2L),
    dest  = c(0L, 2L, 2L, 0L))
res <- assertError(individual_events(events))
check_error(
    res,
    "'time' must be an integer or character vector with non-NA values.")

## Testing animal with only one enter event, keep
events <- data.frame(
    id    = 1L,
    event = 1L,
    time  = 1L,
    node  = 1L,
    dest  = 0L)

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = 1L,
    event = 1L,
    time  = 1L,
    node  = 1L,
    dest  = NA_integer_)

stopifnot(identical(events_obs, events_exp))

## Testing animal with only one exit event, keep
events <- data.frame(
    id    = 1L,
    event = 0L,
    time  = 1L,
    node  = 1L,
    dest  = 0L)

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = 1L,
    event = 0L,
    time  = 1L,
    node  = 1L,
    dest  = NA_integer_)

stopifnot(identical(events_obs, events_exp))

## Testing animal with only one external transfer event, keep
events <- data.frame(
    id    = 1L,
    event = 3L,
    time  = 1L,
    node  = 1L,
    dest  = 2L)

stopifnot(identical(events, as.data.frame(individual_events(events))))

## Testing animal with two enter events, keep first
events <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 1L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = 1L,
    event = 1L,
    time  = 1L,
    node  = 1L,
    dest  = NA_integer_)

stopifnot(identical(events_obs, events_exp))

## Testing animal with two exit events, keep first
events <- data.frame(
    id    = c(1L, 1L),
    event = c(0L, 0L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = 1L,
    event = 0L,
    time  = 1L,
    node  = 1L,
    dest  = NA_integer_)

stopifnot(identical(events_obs, events_exp))

## Testing animal with two enter events and exit, keep path
events <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 1L, 0L),
    time  = c(1L, 2L, 3L),
    node  = c(1L, 2L, 2L),
    dest  = c(0L, 0L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(2L, 3L),
    node  = c(2L, 2L),
    dest  = c(NA_integer_, NA_integer_))

stopifnot(identical(events_obs, events_exp))

## Testing animal with two enter events, a movement and an exit, keep
## path
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 1L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 2L, 1L, 3L),
    dest  = c(0L, 0L, 3L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(1L, 3L, 4L),
    node  = c(1L, 1L, 3L),
    dest  = c(NA_integer_, 3L, NA_integer_))

stopifnot(identical(events_obs, events_exp))

## Testing animal with one enter event and two exit events, keep path
events <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 0L, 0L),
    time  = c(1L, 2L, 3L),
    node  = c(1L, 2L, 1L),
    dest  = c(0L, 0L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(1L, 3L),
    node  = c(1L, 1L),
    dest  = c(NA_integer_, NA_integer_))

stopifnot(identical(events_obs, events_exp))

## Testing animal with one enter event and two exit events, keep path
events <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 0L, 0L),
    time  = c(1L, 2L, 3L),
    node  = c(1L, 1L, 2L),
    dest  = c(0L, 0L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(NA_integer_, NA_integer_))

stopifnot(identical(events_obs, events_exp))

## Testing animal with another event after exit event, exit event
## should be last
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 3L, 0L, 3L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 1L, 2L, 2L),
    dest  = c(0L, 2L, 0L, 3L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(1L, 2L, 3L),
    node  = c(1L, 1L, 2L),
    dest  = c(NA_integer_, 2L, NA_integer_))

stopifnot(identical(events_obs, events_exp))

## Testing animal with another event after exit event,
## no path to exit, don't keep
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(1L, 3L, 0L, 3L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 1L, 3L, 2L),
    dest  = c(0L, 2L, 0L, 3L))

events_obs <- as.data.frame(individual_events(events))

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

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(1L, 3L, 0L),
    time  = c(2L, 3L, 4L),
    node  = c(1L, 1L, 2L),
    dest  = c(NA_integer_, 2L, NA_integer_))

stopifnot(identical(events_obs, events_exp))

pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(individual_events(events))
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Testing animal with another event before enter event, keep path if
## starting on enter event and ending with exit
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(3L, 1L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 2L, 1L, 2L),
    dest  = c(2L, 0L, 2L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(2L, 4L),
    node  = c(2L, 2L),
    dest  = c(NA_integer_, NA_integer_))

stopifnot(identical(events_obs, events_exp))

## Testing animal with no path from enter to exit event, don't keep
events <- data.frame(
    id    = c(1L, 1L, 1L, 1L),
    event = c(3L, 1L, 3L, 0L),
    time  = c(1L, 2L, 3L, 4L),
    node  = c(1L, 2L, 2L, 3L),
    dest  = c(2L, 0L, 1L, 0L))

events_obs <- as.data.frame(individual_events(events))

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

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L, 1L),
    event = c(3L, 3L, 0L),
    time  = c(1L, 2L, 3L),
    node  = c(1L, 2L, 1L),
    dest  = c(2L, 1L, NA_integer_))

stopifnot(identical(events_obs, events_exp))

## Testing animal with no enter or exit event, keep path
events <- data.frame(
    id    = c(1L, 1L),
    event = c(3L, 3L),
    time  = c(1L, 2L),
    node  = c(1L, 2L),
    dest  = c(2L, 3L))

events_obs <- as.data.frame(individual_events(events))

stopifnot(identical(events_obs, events))

## Testing animal with no exit event, keep path
events <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 3L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 2L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 3L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(NA_integer_, 2L))

stopifnot(identical(events_obs, events_exp))

## Testing animal with only enter and exit event, keep path
events <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = c(1L, 1L),
    event = c(1L, 0L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(NA_integer_, NA_integer_))

stopifnot(identical(events_obs, events_obs))

## Testing animal with enter and exit event in wrong order, don't keep
events <- data.frame(
    id    = c(1L, 1L),
    event = c(0L, 1L),
    time  = c(1L, 2L),
    node  = c(1L, 1L),
    dest  = c(0L, 0L))

events_obs <- as.data.frame(individual_events(events))

events_exp <- data.frame(
    id    = integer(0),
    event = integer(0),
    time  = integer(0),
    node  = integer(0),
    dest  = integer(0))

stopifnot(identical(events_obs, events_exp))

## Check converting individual events to u0
events <- data.frame(
    id    = c(1, 1, 1, 1,
              2, 2, 2, 2),
    event = c(1, 3, 3, 0,
              1, 3, 3, 0),
    time  = c(1, 2, 3, 4,
              2, 3, 4, 5),
    node  = c(10, 10, 20, 20,
              10, 10, 20, 20),
    dest  = c(NA, 20, 20, NA,
              NA, 20, 20, NA))

stopifnot(identical(
    u0(individual_events(events), time = 0),
    data.frame(key = c(10, 20),
               node = c(1L, 2L),
               S_1 = c(0L, 0L))))

stopifnot(identical(
    u0(individual_events(events), time = 1),
    data.frame(key = c(10, 20),
               node = c(1L, 2L),
               S_1 = c(1L, 0L))))

stopifnot(identical(
    u0(individual_events(events), time = 2),
    data.frame(key = c(10, 20),
               node = c(1L, 2L),
               S_1 = c(1L, 1L))))

stopifnot(identical(
    u0(individual_events(events), time = 3),
    data.frame(key = c(10, 20),
               node = c(1L, 2L),
               S_1 = c(0L, 2L))))

stopifnot(identical(
    u0(individual_events(events), time = 4),
    data.frame(key = c(10, 20),
               node = c(1L, 2L),
               S_1 = c(0L, 1L))))

stopifnot(identical(
    u0(individual_events(events), time = 5),
    data.frame(key = c(10, 20),
               node = c(1L, 2L),
               S_1 = c(0L, 0L))))

stopifnot(identical(
    u0(individual_events(events), time = 3, age = 2),
    data.frame(key = c(10, 20),
               node = c(1L, 2L),
               S_1 = c(0L, 1L),
               S_2 = c(0L, 1L))))

stopifnot(identical(
    u0(individual_events(events), time = 3, age = 5),
    data.frame(key = c(10, 20),
               node = c(1L, 2L),
               S_1 = c(0L, 2L),
               S_2 = c(0L, 0L))))

stopifnot(identical(
    u0(individual_events(events), time = 3, age = 1),
    data.frame(key = c(10, 20),
               node = c(1L, 2L),
               S_1 = c(0L, 0L),
               S_2 = c(0L, 2L))))

res <- assertError(u0(individual_events(events), time = 4.3))
check_error(
    res,
    "'time' must be an integer or date.")

res <- assertError(u0(individual_events(events),
                      time = c("2021-01-01", "2022-01-01")))
check_error(
    res,
    "'time' must be an integer or date.")

res <- assertError(u0(individual_events(events), time = "2021-01-01"))
check_error(
    res,
    "'time' must be an integer.")

res <- assertError(u0(individual_events(events), time = list()))
check_error(
    res,
    "'time' must be an integer or date.")

res <- assertError(u0(individual_events(events), time = 3, age = -1))
check_error(
    res,
    "'age' must be an integer vector with values > 0.")

events <- data.frame(
    id    = c("individual-1", "individual-1", "individual-1", "individual-1",
              "individual-2", "individual-2", "individual-2", "individual-2"),
    event = c("enter", "extTrans", "extTrans", "exit",
              "enter", "extTrans", "extTrans", "exit"),
    time  = c("2019-02-02", "2020-03-07", "2021-04-14", "2022-05-11",
              "2019-02-02", "2020-03-07", "2021-04-14", "2022-05-11"),
    node  = c("node-1", "node-1", "node-2", "node-2",
              "node-1", "node-1", "node-2", "node-2"),
    dest  = c(NA, "node-2", "node-2", NA,
              NA, "node-2", "node-2", NA))

u0_obs <- u0(individual_events(events))

u0_exp <- data.frame(
    key = c("node-1", "node-2"),
    node = c(1L, 2L),
    S_1 = c(2L, 0L))

stopifnot(identical(u0_obs, u0_exp))

u0_obs <- u0(individual_events(events[rev(seq_len(nrow(events))), ]))

stopifnot(identical(u0_obs, u0_exp))

stopifnot(identical(
    get_individuals(individual_events(events), "2019-02-02"),
    data.frame(
        id = c("individual-1", "individual-2"),
        node = c("node-1", "node-1"),
        age = c(0L, 0L))))

stopifnot(identical(
    get_individuals(individual_events(events), "2019-02-04"),
    data.frame(
        id = c("individual-1", "individual-2"),
        node = c("node-1", "node-1"),
        age = c(2L, 2L))))

stopifnot(identical(
    get_individuals(individual_events(events), "2019-02-01"),
    data.frame(
        id = character(0),
        node = logical(0),
        age = integer(0))))

show_expected <- c(
    "Number of individuals: 2",
    "Number of events: 6")
show_observed <- capture.output(show(individual_events(events)))
stopifnot(identical(show_observed, show_expected))

summary_expected <- c(
    "Number of individuals: 2",
    "Number of events: 6",
    " - Exit: 2",
    " - Enter: 2",
    " - Internal transfer: 0",
    " - External transfer: 2")
summary_observed <- capture.output(summary(individual_events(events)))
stopifnot(identical(summary_observed, summary_expected))

events <- data.frame(
    id    = c(1, 1),
    event = c("extTrans", "exit"),
    time  = c(2, 3),
    node  = c(1, 2),
    dest  = c(2, 0))
res <- assertError(get_individuals(individual_events(events)))
check_error(
    res,
    "All individuals must have an 'enter' event.")

res <- assertError(SimInf:::check_indiv_events_id(3.2))
check_error(
    res,
    "'id' must be an integer or character vector with non-NA values.")

res <- assertError(SimInf:::check_indiv_events_id(NULL))
check_error(
    res,
    "'id' must be an integer or character vector with non-NA values.")
