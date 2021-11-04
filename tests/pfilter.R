## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2021 Stefan Widgren
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

## Create an SIR model object.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5,
             beta = 0.16,
             gamma = 0.077)

## Check that a missing 'time' column in data raises an error.
res <- assertError(SimInf:::pfilter_data(model, data.frame()))
check_error(res, "Missing 'time' column in data.")

## Check that a non-integer value in the 'data$time' column raises an
## error.
res <- assertError(SimInf:::pfilter_data(model, data.frame(time = 1.1)))
check_error(res, "'data$time' must be integer.")

## Check that a NA value in the 'data$time' column raises an error.
res <- assertError(SimInf:::pfilter_data(model, data.frame(time = NA)))
check_error(res, "'data$time' must be integer.")

## Check that data$time[1] < model@tspan[1] raises an error.
res <- assertError(SimInf:::pfilter_data(model, data.frame(time = 0)))
check_error(res, "data$time[1] must be >= tspan[1].")

data <- SimInf:::pfilter_data(model, data.frame(time = 1:3))
stopifnot(identical(
    SimInf:::pfilter_tspan(model, data),
    structure(c(NA, NA, NA, 1, 2, 3), .Dim = 3:2)))

data <- SimInf:::pfilter_data(model, data.frame(time = 2:3))
stopifnot(identical(
    SimInf:::pfilter_tspan(model, data),
    structure(c(1, NA, 2, 3), .Dim = c(2L, 2L))))

## Create an SIR model object where tspan is specified as Dates.
model <- SIR(
    u0 = data.frame(S = 99, I = 1, R = 0),
    tspan = seq(as.Date("2021-01-05"), as.Date("2021-01-09"), by = 1),
    beta = 0.16,
    gamma = 0.077)

## Check that data$time[1] < model@tspan[1] raises an error.
df <- data.frame(time = c("2021-01-04", "2021-01-05"))
res <- assertError(SimInf:::pfilter_data(model, df))
check_error(res, "data$time[1] must be >= tspan[1].")

df <- data.frame(time = c("2021-01-05", "2021-01-06", "2021-01-07"))
data <- SimInf:::pfilter_data(model, df)
stopifnot(identical(
    SimInf:::pfilter_tspan(model, data),
    structure(c(NA, NA, NA, 5, 6, 7), .Dim = 3:2)))

df <- data.frame(time = c("2021-01-06", "2021-01-07"))
data <- SimInf:::pfilter_data(model, df)
stopifnot(identical(
    SimInf:::pfilter_tspan(model, data),
    structure(c(5, NA, 6, 7), .Dim = c(2L, 2L))))

## Check invalid npart
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5,
             beta = 0.16,
             gamma = 0.077)
res <- assertError(pfilter(model = model,
                           npart = 1,
                           data = data.frame(time = 1:3)))
check_error(res, "'npart' must be an integer > 1.")

res <- assertError(pfilter(model = model,
                           npart = c(10, 10),
                           data = data.frame(time = 1:3)))
check_error(res, "'npart' must be an integer > 1.")

## Check the C utility function to split events.
res <- assertError(.Call(SimInf:::SimInf_split_events, 1, 1L))
check_error(res, "'t' must be an integer vector with length >= 1.")

res <- assertError(.Call(SimInf:::SimInf_split_events, 1L, 1))
check_error(res, "'t_end' must be an integer vector with length >= 1.")

stopifnot(identical(.Call(SimInf:::SimInf_split_events, 1L, 1L),
                    structure(c(1L, 1L), .Dim = 1:2)))

stopifnot(identical(.Call(SimInf:::SimInf_split_events, 1L, 1:3),
                    structure(c(1L, 0L, 0L, 1L, 0L, 0L), .Dim = 3:2)))

stopifnot(identical(.Call(SimInf:::SimInf_split_events, 2L, 1:3),
                    structure(c(0L, 1L, 0L, 0L, 1L, 0L), .Dim = 3:2)))

stopifnot(identical(.Call(SimInf:::SimInf_split_events, 3L, 1:3),
                    structure(c(0L, 0L, 1L, 0L, 0L, 1L), .Dim = 3:2)))

stopifnot(identical(.Call(SimInf:::SimInf_split_events, 4L, 1:3),
                    structure(c(0L, 0L, 0L, 0L, 0L, 0L), .Dim = 3:2)))

## Split events
events <- data.frame(
    event      = rep("extTrans", 6),
    time       = c(1, 1, 2, 2, 3, 3),
    node       = c(3, 3, 1, 4, 3, 4),
    dest       = c(4, 2, 3, 3, 2, 2),
    n          = c(9, 2, 8, 3, 5, 4),
    proportion = c(0, 0, 0, 0, 0, 0),
    select     = c(4, 4, 4, 4, 4, 4),
    shift      = c(0, 0, 0, 0, 0, 0))

model <- SIR(u0 = data.frame(S = c(10, 15, 20, 25),
                             I = c(0, 0, 0, 0),
                             R = c(0, 0, 0, 0)),
             tspan = 0:3,
             beta = 0.16,
             gamma = 0.077,
             events = events)

stopifnot(identical(
    SimInf:::pfilter_events(model@events, 2:3),
    list(new("SimInf_events",
             E = new("dgCMatrix",
                     i = c(0L, 1L, 2L, 0L, 1L, 2L),
                     p = c(0L, 1L, 2L, 3L, 6L),
                     Dim = 3:4,
                     Dimnames = list(c("S", "I", "R"), c("1", "2", "3", "4")),
                     x = c(1, 1, 1, 1, 1, 1),
                     factors = list()),
             N = structure(integer(0),
                           .Dim = c(0L, 0L)),
             event = c(3L, 3L, 3L, 3L),
             time = c(1L, 1L, 2L, 2L),
             node = c(3L, 3L, 1L, 4L),
             dest = c(4L, 2L, 3L, 3L),
             n = c(9L, 2L, 8L, 3L),
             proportion = c(0, 0, 0, 0),
             select = c(4L, 4L, 4L, 4L),
             shift = c(0L, 0L, 0L, 0L)),
         new("SimInf_events",
             E = new("dgCMatrix",
                     i = c(0L, 1L, 2L, 0L, 1L, 2L),
                     p = c(0L, 1L, 2L, 3L, 6L),
                     Dim = 3:4,
                     Dimnames = list(c("S", "I", "R"), c("1", "2", "3", "4")),
                     x = c(1, 1, 1, 1, 1, 1),
                     factors = list()),
             N = structure(integer(0),
                           .Dim = c(0L, 0L)),
             event = c(3L, 3L),
             time = c(3L, 3L),
             node = 3:4,
             dest = c(2L, 2L),
             n = 5:4,
             proportion = c(0, 0),
             select = c(4L, 4L),
             shift = c(0L, 0L)))))

## Check that an error is raised if a weight is invalid.
res <- assertError(.Call(SimInf:::SimInf_systematic_resampling, NaN))
check_error(res, "Invalid weight detected (non-finite or < 0.0).")
res <- assertError(.Call(SimInf:::SimInf_systematic_resampling, -0.1))
check_error(res, "Invalid weight detected (non-finite or < 0.0).")

## Check that sum of weights >= 0.
res <- assertError(.Call(SimInf:::SimInf_systematic_resampling, 0))
check_error(res, "Non-positive sum of weights detected.")

## Expect all particles if the weights are equal.
w <- rep(0.1, 10)
stopifnot(identical(
    seq_len(length(w)),
    .Call(SimInf:::SimInf_systematic_resampling, w)))

## Expect function for the observation process
obs_fn <- function(model, data) {0}

stopifnot(identical(
    SimInf:::pfilter_obs_process(
                 SIR(u0 = data.frame(S = 99, I = 1, R = 0),
                     tspan = seq(1, 21, by = 3),
                     beta = 0.16,
                     gamma = 0.077),
                 obs_fn,
                 list("1" = data.frame(time = 1, Iobs = 1),
                      "4" = data.frame(time = 4, Iobs = 2),
                      "7" = data.frame(time = 7, Iobs = 2),
                      "10" = data.frame(time = 10, Iobs = 3),
                      "13" = data.frame(time = 13, Iobs = 3),
                      "16" = data.frame(time = 16, Iobs = 3),
                      "19" = data.frame(time = 19, Iobs = 3)),
                 5),
    obs_fn))

## Raise an error if the observation process is not a function when a
## model contains multiple nodes.
res <- assertError(
    SimInf:::pfilter_obs_process(
                 SIR(u0 = data.frame(S = c(99, 99), I = 1, R = 0),
                     tspan = seq(1, 21, by = 3),
                     beta = 0.16,
                     gamma = 0.077),
                 Iobs ~ poisson(I + 1e-6),
                 list("1" = data.frame(time = 1, Iobs = 1),
                      "4" = data.frame(time = 4, Iobs = 2),
                      "7" = data.frame(time = 7, Iobs = 2),
                      "10" = data.frame(time = 10, Iobs = 3),
                      "13" = data.frame(time = 13, Iobs = 3),
                      "16" = data.frame(time = 16, Iobs = 3),
                      "19" = data.frame(time = 19, Iobs = 3)),
                 5))

check_error(res,
            paste("The observation process must be a function",
                  "for a model with multiple nodes."))

## Raise an error if the observation process is not a function when
## data contains multiple rows for a time-point.
res <- assertError(
    SimInf:::pfilter_obs_process(
                 SIR(u0 = data.frame(S = 99, I = 1, R = 0),
                     tspan = seq(1, 21, by = 3),
                     beta = 0.16,
                     gamma = 0.077),
                 Iobs ~ poisson(I + 1e-6),
                 list("1" = data.frame(time = c(1, 1), Iobs = c(1, 2)),
                      "7" = data.frame(time = 7, Iobs = 2),
                      "10" = data.frame(time = 10, Iobs = 3),
                      "13" = data.frame(time = 13, Iobs = 3),
                      "16" = data.frame(time = 16, Iobs = 3),
                      "19" = data.frame(time = 19, Iobs = 3)),
                 5))

check_error(res,
            paste("The observation process must be a function",
                  "when data contains multiple rows for a time-point."))

## Raise an error if the observation process is not a function or a
## formula.
res <- assertError(
    SimInf:::pfilter_obs_process(
                 SIR(u0 = data.frame(S = 99, I = 1, R = 0),
                     tspan = seq(1, 21, by = 3),
                     beta = 0.16,
                     gamma = 0.077),
                 "obs_process",
                 list("1" = data.frame(time = 1, Iobs = 1),
                      "4" = data.frame(time = 4, Iobs = 2),
                      "7" = data.frame(time = 7, Iobs = 2),
                      "10" = data.frame(time = 10, Iobs = 3),
                      "13" = data.frame(time = 13, Iobs = 3),
                      "16" = data.frame(time = 16, Iobs = 3),
                      "19" = data.frame(time = 19, Iobs = 3)),
                 5))

check_error(res, "'obs_process' must be either a formula or a function.")

## Raise an error if the lhs does not match a column in data.
res <- assertError(
    SimInf:::pfilter_obs_process(
                 SIR(u0 = data.frame(S = 99, I = 1, R = 0),
                     tspan = seq(1, 21, by = 3),
                     beta = 0.16,
                     gamma = 0.077),
                 Robs ~ poisson(I + 1e-6),
                 list("1" = data.frame(time = 1, Iobs = 1),
                      "4" = data.frame(time = 4, Iobs = 2),
                      "7" = data.frame(time = 7, Iobs = 2),
                      "10" = data.frame(time = 10, Iobs = 3),
                      "13" = data.frame(time = 13, Iobs = 3),
                      "16" = data.frame(time = 16, Iobs = 3),
                      "19" = data.frame(time = 19, Iobs = 3)),
                 5))

check_error(res, "Unable to match the parameter on the lhs to a column in 'data'.")
