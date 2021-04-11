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

##' Split tspan into intervals.
##' @noRd
pfilter_tspan <- function(model, data) {
    if (!is.data.frame(data))
        data <- as.data.frame(data)

    if (!("time" %in% names(data)))
        stop("Missing 'time' column in data.", call. = FALSE)

    if (!is.null(names(model@tspan))) {
        data$time <- julian(x = as.Date(data$time),
                            origin = as.Date(names(model@tspan)[1]))
        data$time <- as.numeric(data$time + model@tspan[1])
    }

    check_integer_arg(data$time)

    if (any(length(data$time) < 1,
            any(diff(data$time) <= 0),
            any(is.na(data$time)))) {
        stop("'time' column in data must be an increasing vector.",
             call. = FALSE)
    }

    if (data$time[1] < model@tspan[1])
        stop("data$time[1] must be >= tspan[1].", call. = TRUE)

    do.call("rbind", lapply(seq_len(length(data$time)), function(i) {
        if (i == 1) {
            if (model@tspan[1] < data$time[1])
                return(as.numeric(c(model@tspan[1], data$time[1])))
            return(c(NA_real_, data$time[i]))
        }

        if (diff(as.integer(data$time[c(i - 1L, i)])) > 1)
            return(as.numeric(data$time[c(i - 1L, i)]))
        c(NA_real_, data$time[i])
    }))
}
