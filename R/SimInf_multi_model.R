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
##' @slot multi_model FIXME.
##' @slot events FIXME.
##' @slot data A \code{list} with \code{data.frame} items holding the
##'     time series data.
##' @slot n_models Integer with the number of models.
##' @export
setClass(
    "SimInf_multi_model",
    slots = c(model       = "SimInf_model",
              multi_model = "SimInf_model",
              events      = "list",
              data        = "list",
              n_models    = "integer")
)

##' Create a multi-model object
##'
##' @param model The \code{SimInf_model} object to estimate parameters
##'     in.
##' @param multi_model FIXME.
##' @param data A \code{data.frame} holding the time series data.
##' @value FIXME
##' @export
multi_model <- function(model, multi_model, data) {
    if (any(isFALSE(identical(dim(model@U_sparse), c(0L, 0L))),
            isFALSE(identical(dim(model@V_sparse), c(0L, 0L))),
            isFALSE(identical(dim(multi_model@U_sparse), c(0L, 0L))),
            isFALSE(identical(dim(multi_model@V_sparse), c(0L, 0L))))) {
        stop("Cannot create a multi model object with a sparse result matrix.",
             call. = FALSE)
    }

    if (any(n_nodes(multi_model) <= n_nodes(model),
            n_nodes(multi_model) %% n_nodes(model))) {
        stop("Invalid number of nodes in the multi_model object.",
             call. = FALSE)
    }

    data <- pfilter_data(multi_model, data)
    tspan <- pfilter_tspan(model, data)
    multi_model@tspan <- tspan[, 2]
    events <- pfilter_events(multi_model@events, tspan[, 2])
    n_models <- as.integer(n_nodes(multi_model) / n_nodes(model))

    methods::new("SimInf_multi_model",
                 model = model,
                 multi_model = multi_model,
                 events = events,
                 data = data,
                 n_models = n_models)
}
