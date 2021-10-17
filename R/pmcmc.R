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
##' @slot npart n integer with the number of particles (> 1) to use in
##'     the bootstrap particle filter.
##' @slot obs_process A \code{formula} or \code{function} determining
##'     the observation process.
##' @slot data A \code{data.frame} holding the time series data for
##'     the observation process.
##' @slot chain FIXME
##' @slot pf FIXME
##' @export
setClass(
    "SimInf_pmcmc",
    slots = c(model       = "SimInf_model",
              priors      = "data.frame",
              target      = "character",
              pars        = "integer",
              npart       = "integer",
              obs_process = "ANY",
              data        = "data.frame",
              chain       = "matrix",
              pf          = "list",
              adaptmix    = "numeric")
)

##' Check if a SimInf_pmcmc object is valid
##'
##' @param object The SimInf_pmcmc object to check.
##' @noRd
valid_SimInf_pmcmc_object <- function(object) {
    errors <- character(0)

    if (length(object@adaptmix) != 1L ||
        object@adaptmix <= 0 ||
        object@adaptmix >= 1) {
        errors <- c(errors, "'adaptmix' must be a value > 0 and < 1")
    }

    if (!identical(object@target, "gdata") &&
        !identical(object@target, "ldata")) {
        errors <- c(errors, "'target' must be 'gdata' or 'ldata'")
    }

    if (length(errors))
        return(errors)
    TRUE
}

## Assign the validity method for the SimInf_pmcmc class.
setValidity("SimInf_pmcmc", valid_SimInf_pmcmc_object)

summary_chain <- function(chain) {
    qq <- do.call("rbind", apply(chain, 2, function(x) {
        cbind(t(quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))),
              Mean = mean(x),
              SD = sqrt(var(x, na.rm = TRUE)))
    }, simplify = FALSE))
    rownames(qq) <- colnames(chain)
    print.table(qq, digits = 3)
}

##' Brief summary of a \code{SimInf_pmcmc} object
##'
##' @param object The \code{SimInf_pmcmc} object.
##' @return \code{invisible(object)}.
##' @export
##' @importFrom methods show
setMethod(
    "show",
    signature(object = "SimInf_pmcmc"),
    function(object) {
        cat("Particle Markov chain Monte Carlo\n")
        cat("---------------------------------\n")
        cat(sprintf("Number of iterations: %i\n", length(object)))
        cat(sprintf("Number of particles: %i\n", object@npart))
        cat(sprintf("Mixing proportion for adaptive proposal: %.2f\n",
                    object@adaptmix))

        if (length(object) > 0) {
            cat(sprintf("Acceptance ratio: %.3f\n",
                        mean(object@chain[, "accept"])))

            print_title(
                "Quantiles, mean and standard deviation for each variable")
            summary_chain(object@chain[, -(1:4)])
        }

        invisible(object)
    }
)

##' Particle Markov chain Monte Carlo (PMCMC) algorithm
##'
##' @param model The model to simulate data from.
##' @template obs_process-param
##' @template data-param
##' @template priors-param
##' @template npart-param
##' @param niter An integer specifying the number of iterations to run
##'     the PMCMC.
##' @param adaptmix Mixing proportion for adaptive proposal.
##' @template verbose-param
##' @references
##'
##' \Andrieu2010
##' @export
setGeneric(
    "pmcmc",
    signature = "model",
    function(model, obs_process, data, priors, npart, niter,
             theta = NULL, adaptmix = 0.05,
             verbose = getOption("verbose", FALSE)) {
        standardGeneric("pmcmc")
    }
)

##' @rdname pmcmc
##' @export
setMethod(
    "pmcmc",
    signature(model = "SimInf_model"),
    function(model, obs_process, data, priors, npart, niter, theta,
             adaptmix, verbose) {
        check_integer_arg(npart)
        npart <- as.integer(npart)
        if (length(npart) != 1L || npart <= 1L)
            stop("'npart' must be an integer > 1.", call. = FALSE)

        adaptmix <- as.numeric(adaptmix)
        if (length(adaptmix) != 1L || adaptmix <= 0 || adaptmix >= 1)
            stop("'adaptmix' must be a value > 0 and < 1", call. = FALSE)

        ## Match the 'priors' to parameters in 'ldata' or 'gdata'.
        priors <- parse_priors(priors)
        pars <- match_priors(model, priors)

        object <- new("SimInf_pmcmc", model = model, priors = priors,
                      target = pars$target, pars = pars$pars,
                      obs_process = obs_process, data = data,
                      npart = npart, adaptmix = adaptmix)

        if (!is.null(theta)) {
            check_integer_arg(niter)
            niter <- as.integer(niter)
            if (length(niter) != 1L || niter <= 0L)
                stop("'niter' must be an integer > 0.", call. = FALSE)

            if (!all(is.atomic(theta),
                     is.numeric(theta),
                     all(priors$parameter %in% names(theta)))) {
                stop("'theta' must be a vector with initial ",
                     "values for the parameters.",
                     call. = FALSE)
            }
            theta <- theta[priors$parameter]

            npars <- length(object@pars)
            object@chain <- setup_chain(object, 1L)
            object@pf <- setup_pf(object, 1L)

            if (object@target == "gdata") {
                for (i in seq_len(npars)) {
                    object@model@gdata[object@pars[i]] <- theta[i]
                }
            } else {
                for (i in seq_len(npars)) {
                    object@model@ldata[object@pars[i], ] <- theta[i]
                }
            }

            object@pf[[1]] <- pfilter(object@ model,
                                      obs_process = object@obs_process,
                                      object@data,
                                      npart = object@npart)

            loglik <- object@pf[[1]]@loglik
            logprior <- dpriors(theta, object@priors)
            logpost <- loglik + logprior
            accept <- 0
            object@chain[1, ] <- c(logpost, loglik, logprior, accept, theta)

            niter <- niter - 1L
            if (niter == 0)
                return(object)
        }

        continue(object, niter = niter, verbose = verbose)
    }
)

is_empty_chain <- function(object) {
    isTRUE(length(object) == 0L)
}

setup_chain <- function(object, niter) {
    m <- matrix(NA_real_,
                nrow = niter,
                ncol = 4L + length(object@pars),
                dimnames = list(NULL, c("logpost", "loglik",
                                        "logprior", "accept",
                                        object@priors$parameter)))

    if (is_empty_chain(object))
        return(m)
    rbind(object@chain, m)
}

setup_pf <- function(object, niter) {
    c(object@pf, lapply(seq_len(niter), function(i) NULL))
}

pmcmc_progress <- function(object, i, verbose) {
    if (isTRUE(verbose) && isTRUE(i %% 100 == 0)) {
        print_title(sprintf(
            "PMCMC iteration: %i of %i. Acceptance ratio: %.3f",
            i, length(object),
            mean(object@chain[seq_len(i), "accept"])))
        summary_chain(object@chain[seq_len(i), vars_i])
    }

    invisible(NULL)
}

##' Length of the MCMC chain
##'
##' @param x The \code{SimInf_pmcmc} object determine the length of
##'     the MCMC chain for.
##' @export
setMethod(
    "length",
    signature(x = "SimInf_pmcmc"),
    function(x) {
        nrow(x@chain)
    }
)

##' @rdname continue
##' @importFrom mvtnorm rmvnorm
##' @export
setMethod(
    "continue",
    signature(object = "SimInf_pmcmc"),
    function(object, niter, ...,
             verbose = getOption("verbose", FALSE)) {
        check_integer_arg(niter)
        niter <- as.integer(niter)
        if (length(niter) != 1L || niter <= 0L)
            stop("'niter' must be an integer > 0.", call. = FALSE)

        iterations <- length(object) + seq_len(niter)
        npars <- length(object@pars)
        object@chain <- setup_chain(object, niter)
        object@pf <- setup_pf(object, niter)

        if (iterations[1] == 1L) {
            iterations <- iterations[-1]

            theta <- rpriors(object@priors)
            if (object@target == "gdata") {
                for (i in seq_len(npars)) {
                    object@model@gdata[object@pars[i]] <- theta[i]
                }
            } else {
                for (i in seq_len(npars)) {
                    object@model@ldata[object@pars[i], ] <- theta[i]
                }
            }

            object@pf[[1]] <- pfilter(object@ model,
                                      obs_process = object@obs_process,
                                      object@data,
                                      npart = object@npart)

            pf <- object@pf[[1]]
            loglik <- pf@loglik
            logprior <- dpriors(theta, object@priors)
            logpost <- loglik + logprior
            accept <- 0
            object@chain[1, ] <- c(logpost, loglik, logprior, accept, theta)
        } else {
            ## Continue from the last iteration in the chain.
            i <- iterations[1] - 1
            pf <- object@pf[[i]]
            logpost <- object@chain[i, "logpost"]
            loglik <- object@chain[i, "loglik"]
            logprior <- object@chain[i, "logprior"]
            theta <- object@chain[i, -(1:4)]
        }

        for (i in iterations) {
            ## Proposal
            if (runif(1) < object@adaptmix || i <= 2 * npars) {
                sigma <- diag(0.1^2 / npars, npars)
            } else if (npars == 1) {
                sigma <- matrix(2.38^2 * var(object@chain[seq_len(i - 1), -c(1:4)]))
            } else {
                sigma <- 2.38^2 / npars * cov(object@chain[seq_len(i - 1), -c(1:4)])
            }
            theta_prop <- rmvnorm(n = 1, mean = theta, sigma = sigma)[1, ]
            logprior_prop <- dpriors(theta_prop, object@priors)

            if (is.finite(logprior_prop)) {
                if (object@target == "gdata") {
                    for (j in seq_len(npars)) {
                        object@model@gdata[object@pars[j]] <- theta_prop[j]
                    }
                } else {
                    for (j in seq_len(npars)) {
                        object@model@ldata[object@pars[j], ] <- theta_prop[j]
                    }
                }

                pf_prop <- pfilter(object@model,
                                   obs_process = object@obs_process,
                                   object@data,
                                   npart = object@npart)
                loglik_prop <- pf_prop@loglik

                accept <- 0
                alpha <- exp(loglik_prop + logprior_prop - loglik - logprior)
                if (is.finite(alpha) && runif(1) < alpha) {
                    loglik <- loglik_prop
                    logprior <- logprior_prop
                    logpost <- loglik + logprior
                    theta <- theta_prop
                    pf <- pf_prop
                    accept <- 1
                }
            }

            ## Save current value of chain.
            object@chain[i, ] <- c(logpost, loglik, logprior, accept, theta)
            object@pf[[i]] <- pf

            ## Report progress.
            pmcmc_progress(object, i, verbose)
        }

        object
    }
)
