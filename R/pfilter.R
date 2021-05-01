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

##' Class \code{"SimInf_pfilter"}
##'
##' @slot model The \code{SimInf_model} object to simulate data from.
##' @slot npart An integer with the number of particles that was used
##'     at each timestep.
##' @slot tspan A two-column matrix where each row specifies the tspan
##'     to use when running the model from time[i] to time[i+1]. The
##'     first column is \code{NA} if the interval is one time-unit.
##' @slot loglik The estimated log likelihood.
##' @slot ess A numeric vector with the effective sample size (ESS).
##'     The effective sample size is computed as
##'     \deqn{\left(\sum_{i=1}^N\!(w_{t}^{(i)})^2\right)^{-1},}{1 /
##'     (sum(w_it^2)),} where \eqn{w_{t}^{(i)}}{w_it} is the
##'     normalized weight of particle \eqn{i} at time \eqn{t}.
##' @export
setClass(
    "SimInf_pfilter",
    slots = c(model  = "SimInf_model",
              npart  = "integer",
              tspan  = "matrix",
              loglik = "numeric",
              ess    = "numeric")
)

##' Print summary of a \code{SimInf_pfilter} object
##'
##' @param object The \code{SimInf_pfilter} object.
##' @return \code{invisible(object)}.
##' @export
##' @importFrom methods show
setMethod(
    "show",
    signature(object = "SimInf_pfilter"),
    function(object) {
        cat(sprintf("Number of particles: %i\n", object@npart))
        cat(sprintf("Log-likelihood: %f\n", object@loglik))

        invisible(object)
    }
)

##' Print summary of a \code{SimInf_pfilter} object
##'
##' @param object The \code{SimInf_pfilter} object.
##' @param ... Unused additional arguments.
##' @return \code{invisible(NULL)}.
##' @export
setMethod(
    "summary",
    signature(object = "SimInf_pfilter"),
    function(object, ...) {
        cat(sprintf("Number of particles: %i\n", object@npart))
        cat(sprintf("Log-likelihood: %f\n", object@loglik))

        invisible(NULL)
    }
)

##' Validate the number of particles.
##' @noRd
pfilter_npart <- function(npart) {
    check_integer_arg(npart)
    npart <- as.integer(npart)
    if (length(npart) != 1L || npart <= 1L)
        stop("'npart' must be an integer > 1.", call. = FALSE)
    npart
}

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

pfilter_obs_process <- function(model, obs_process) {
    if (is.function(obs_process))
        return(match.fun(obs_process))

    if (n_nodes(model) > 1) {
        stop("The observation process must be a function ",
             "for a model with multiple nodes.",
             call. = FALSE)
    }

    if (!is(obs_process, "formula")) {
        stop("'obs_process' must be either a formula or a function.",
             call. = FALSE)
    }

    stop("Not implemented", call. = FALSE)
}

##' Bootstrap particle filter
##'
##' Systematic resampling is performed at each observation.
##'
##' @param model The \code{SimInf_model} object to simulate data from.
##' @param obs_process A \code{formula} or \code{function} determining
##'     the observation process.
##' @param data A \code{data.frame} holding the time series data.
##' @param npart An integer with the number of particles (> 1) to use
##'     at each timestep.
##' @return A \code{SimInf_pfilter} object.
##' @references
##'
##' \Gordon1993
##' @export
setGeneric(
    "pfilter",
    signature = "model",
    function(model, obs_process, data, npart) {
        standardGeneric("pfilter")
    }
)

##' @rdname pfilter
##' @export
setMethod(
    "pfilter",
    signature(model = "SimInf_model"),
    function(model, obs_process, data, npart) {
        npart <- pfilter_npart(npart)
        tspan <- pfilter_tspan(model, data)
        obs_process <- pfilter_obs_process(model, obs_process)

        if (n_nodes(model) == 1)
            return(pfilter_single_node(model, obs_process, data, npart, tspan))
        pfilter_multiple_nodes(model, obs_process, data, npart, tspan)
    }
)
