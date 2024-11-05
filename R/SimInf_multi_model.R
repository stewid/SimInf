## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2024 Stefan Widgren
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

##' Class \code{"SimInf_multi_model"}
##'
##' @slot model The \code{SimInf_model} object to estimate parameters
##'     in.
##' @slot tspan FIXME.
##' @slot events FIXME.
##' @slot data A \code{list} with \code{data.frame} items holding the
##'     time series data.
##' @slot n_models Integer with the number of models.
##' @export
setClass(
    "SimInf_multi_model",
    slots = c(model       = "SimInf_model",
              tspan       = "matrix",
              events      = "list",
              data        = "list",
              n_models    = "integer")
)

##' Create a multi-model object
##'
##' @param model The \code{SimInf_model} object to estimate parameters
##'     in.
##' @param n_models FIXME.
##' @param data A \code{data.frame} holding the time series data.
##' @return FIXME
##' @export
multi_model <- function(model, n_models, data) {
    if (any(isFALSE(identical(dim(model@U_sparse), c(0L, 0L))),
            isFALSE(identical(dim(model@V_sparse), c(0L, 0L))))) {
        stop("Cannot create a multi model object with a sparse result matrix.",
             call. = FALSE)
    }

    data <- pfilter_data(model, data)
    tspan <- pfilter_tspan(model, data)
    model@tspan <- tspan[, 2]
    events <- pfilter_events(model@events, tspan[, 2])
    n_models <- as.integer(n_models)

    methods::new("SimInf_multi_model",
                 model = model,
                 tspan = tspan,
                 events = events,
                 data = data,
                 n_models = n_models)
}
