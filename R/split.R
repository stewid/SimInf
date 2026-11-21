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

##' Split observed data by observation time.
##'
##' Validate the observed data and split it into one data.frame per
##' unique observation time, ordered by increasing time. Rows that
##' share an observation time are kept together in the same
##' data.frame. If the model's tspan has names (i.e. is expressed in
##' dates), the time column in data is converted from dates to
##' numeric time points, relative to the first day of tspan.
##'
##' @param model The \code{SimInf_model} object that determines the
##'     observation time points.
##' @param data A \code{data.frame} with the observed data. It must
##'     contain a \code{time} column with the observation time
##'     points, coerced to dates when tspan is expressed in dates.
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
