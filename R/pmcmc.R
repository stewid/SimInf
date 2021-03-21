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

##' Class \code{"SimInf_pmcmc"}
##'
##' @slot model The \code{SimInf_model} object to estimate parameters
##'     in.
##' @template priors-slot
##' @slot target Character vector (\code{gdata} or \code{ldata}) that
##'     determines if the \code{pmcmc} method estimates parameters in
##'     \code{model@@gdata} or in \code{model@@ldata}.
##' @slot pars Index to the parameters in \code{target}.
##' @export
setClass(
    "SimInf_pmcmc",
    slots = c(model   = "SimInf_model",
              priors  = "data.frame",
              target  = "character",
              pars    = "integer")
)

##' Particle Markov chain Monte Carlo (PMCMC) algorithm
##'
##' @param model The model to simulate data from.
##' @param data A data frame holding the time series data.
##' @template priors-param
##' @param npart An integer specifying the number of particles for the
##'     bootstrap particle filter.
##' @param niter An integer specifying the number of iterations to run
##'     the PMCMC.
##' @param fn FIXME.
##' @param ... Further arguments to be passed to \code{fn}.
##' @template verbose-param
##' @references
##'
##' \Andrieu2010
##' @export
pmcmc <- function(model, data, priors, npart, niter, fn, ...,
                  verbose = getOption("verbose", FALSE)) {
    check_model_argument(model)

    check_integer_arg(npart)
    npart <- as.integer(npart)
    if (length(npart) != 1L || npart <= 1L)
        stop("'npart' must be an integer > 1.", call. = FALSE)

    check_integer_arg(niter)
    niter <- as.integer(niter)
    if (length(niter) != 1L || niter <= 0L)
        stop("'niter' must be an integer > 0.", call. = FALSE)

    ## Match the 'priors' to parameters in 'ldata' or 'gdata'.
    priors <- parse_priors(priors)
    pars <- match_priors(model, priors)

    object <- new("SimInf_pmcmc", model = model, priors = priors,
                  target = pars$target, pars = pars$pars)

    continue(object, niter = niter, verbose = verbose, ...)
}
