## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2026 Stefan Widgren
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

##' Split observed data by observation time
##'
##' Validate the observed data and split it into one data.frame per
##' unique observation time, ordered by increasing time. Rows that
##' share an observation time are kept together in the same
##' data.frame. If the model's tspan has names (i.e. is expressed in
##' dates), the time column in data is converted from dates to numeric
##' time points, relative to the first day of tspan.
##'
##' @param model The \code{SimInf_model} object that determines the
##'     observation time points.
##' @param data A \code{data.frame} with the observed data. It must
##'     contain a \code{time} column with the observation time points,
##'     coerced to dates when tspan is expressed in dates.
##' @return A named \code{list} with one \code{data.frame} per unique
##'     observation time, ordered by increasing time and with the
##'     observation time as the name of each list item.
##' @noRd
split_data <- function(model, data) {
    if (!is.data.frame(data))
        stop("'data' must be a data.frame.", call. = FALSE)

    if (!("time" %in% names(data)))
        stop("Missing 'time' column in data.", call. = FALSE)

    if (!is.null(names(model@tspan))) {
        data$time <- julian(x = as.Date(data$time),
                            origin = as.Date(names(model@tspan)[1]))
        data$time <- as.numeric(data$time + model@tspan[1])
    }

    check_integer_arg(data$time)
    data <- data[order(data$time), seq_len(ncol(data)), drop = FALSE]

    if (any(length(data$time) < 1,
            anyNA(data$time),
            any(diff(unique(data$time)) <= 0))) {
        stop("'time' column in data must be an increasing vector.",
             call. = FALSE)
    }

    if (data$time[1] < model@tspan[1])
        stop("data$time[1] must be >= tspan[1].", call. = FALSE)

    lapply(split(data, data$time), function(x) {
        rownames(x) <- NULL
        x
    })
}

##' Split the observation times into intervals for simulation
##'
##' Determine the 'tspan' to use for each call to \code{run} when
##' propagating particles from one observation time to the next.
##' Consecutive observation times that are more than one time-unit
##' apart define an interval by both its start and end time.  When the
##' gap is exactly one time-unit, only the end time is needed to
##' advance the model and the start time is set to \code{NA}. The
##' first interval starts at the first time point in 'tspan', unless
##' the data already start there, in which case it is a one time-unit
##' interval specified by its end time only.
##'
##' @param model The \code{SimInf_model} object with the 'tspan' that
##'     determines the first time point of the simulation.
##' @param data A named list with one data.frame per observation time,
##'     i.e. the output from \code{split_data}.
##' @return A two-column numeric matrix with one row per observation
##'     time. The second column holds the end time of each interval;
##'     the first column holds the start time, or \code{NA} when the
##'     interval spans a single time-unit and the end time alone
##'     determines the run.
##' @noRd
split_tspan <- function(model, data) {
    time <- sapply(data, function(x) {
        x$time[1]
    })

    do.call("rbind", lapply(seq_along(time), function(i) {
        if (i == 1) {
            if (model@tspan[1] < time[1])
                return(as.numeric(c(model@tspan[1], time[1])))
            return(as.numeric(c(NA_real_, time[i])))
        }

        if (diff(as.integer(time[c(i - 1L, i)])) > 1)
            return(as.numeric(time[c(i - 1L, i)]))
        as.numeric(c(NA_real_, time[i]))
    }))
}
