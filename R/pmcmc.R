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
##' @slot chain A matrix where each row contains \code{logPost},
##'     \code{logLik}, \code{logPrior}, \code{accept}, and the
##'     \code{parameters} for each iteration.
##' @slot pf List with the filtered trajectory from each iteration.
##' @slot adaptmix Mixing proportion for adaptive proposal.
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
        errors <- c(errors, "'adaptmix' must be a value >= 0 and <= 1")
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

            ## Skip first four columns in chain.
            j <- seq(from = 5, by = 1, length.out = length(object@pars))
            summary_chain(object@chain[, j, drop = FALSE])
        }

        invisible(object)
    }
)

##' Detailed summary of a \code{SimInf_pmcmc} object
##'
##' @param object The \code{SimInf_pmcmc} object
##' @param ... Not used.
##' @return None (invisible 'NULL').
##' @export
setMethod(
    "summary",
    signature(object = "SimInf_pmcmc"),
    function(object, ...) {
        cat("Particle Markov chain Monte Carlo\n")
        cat("---------------------------------\n")
        cat(sprintf("Number of iterations: %i\n", length(object)))
        cat(sprintf("Number of particles: %i\n", object@npart))
        if (length(object) > 0) {
            cat(sprintf("Acceptance ratio: %.3f\n",
                        mean(object@chain[, "accept"])))
        }

        ## The model name
        cat(sprintf("Model: %s\n", as.character(class(object@model))))

        ## Nodes
        cat(sprintf("Number of nodes: %i\n", n_nodes(object@model)))

        summary_transitions(object@model)

        if (length(object) > 0) {
            print_title(
                "Quantiles, mean and standard deviation for each variable")

            ## Skip the first four columns in chain: 'logPost',
            ## 'logLik', 'logPrior', and 'accept'.
            j <- seq(from = 5, by = 1, length.out = length(object@pars))
            summary_chain(object@chain[, j, drop = FALSE])
        }

        invisible(NULL)
    }
)

##' Particle Markov chain Monte Carlo (PMCMC) algorithm
##'
##' @param model The model to simulate data from.
##' @template obs_process-param
##' @template data-param
##' @template priors-param
##' @template npart-param
##' @template niter-param
##' @param theta A named vector of initial values for the parameters
##'     of the model.  Default is \code{NULL}, and then these are
##'     sampled from the prior distribution(s).
##' @param adaptmix Mixing proportion for adaptive proposal.  Must be
##'     a value between zero and one.
##' @template verbose-param
##' @references
##'
##' \Andrieu2010
##'
##' \Roberts2009
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
        if (any(length(npart) != 1L,
                any(npart <= 1L))) {
            stop("'npart' must be an integer > 1.", call. = FALSE)
        }

        adaptmix <- as.numeric(adaptmix)
        if (any(length(adaptmix) != 1L,
                any(adaptmix <= 0),
                any(adaptmix >= 1))) {
            stop("'adaptmix' must be a value > 0 and < 1.", call. = FALSE)
        }

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
            if (any(length(niter) != 1L,
                    any(niter <= 0L))) {
                stop("'niter' must be an integer > 0.", call. = FALSE)
            }

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

            methods::slot(object@model, object@target) <-
                set_proposal(object, theta)
            object@pf[[1]] <- pfilter(object@model,
                                      obs_process = object@obs_process,
                                      object@data,
                                      npart = object@npart)

            logLik <- object@pf[[1]]@loglik
            logPrior <- dpriors(theta, object@priors)
            logPost <- logLik + logPrior
            accept <- 0
            object@chain[1, ] <- c(logPost, logLik, logPrior, accept, theta)

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
                dimnames = list(NULL, c("logPost", "logLik",
                                        "logPrior", "accept",
                                        object@priors$parameter)))

    if (is_empty_chain(object))
        return(m)
    rbind(object@chain, m)
}

setup_pf <- function(object, niter) {
    c(object@pf, vector("list", length = niter))
}

set_proposal <- function(object, theta) {
    if (object@target == "gdata") {
        for (i in seq_len(length(object@pars))) {
            object@model@gdata[object@pars[i]] <- theta[i]
        }
    } else {
        for (i in seq_len(length(object@pars))) {
            object@model@ldata[object@pars[i], ] <- theta[i]
        }
    }

    methods::slot(object@model, object@target)
}

pmcmc_progress <- function(object, i, verbose) {
    if (isTRUE(verbose) && isTRUE(i %% 100 == 0)) {
        print_title(sprintf(
            "PMCMC iteration: %i of %i. Acceptance ratio: %.3f",
            i, length(object),
            mean(object@chain[seq_len(i), "accept"])))

        ## Skip columns logLik, logPrior and accept in the chain.
        j <- c(1, seq(from = 5, by = 1, length.out = length(object@pars)))
        summary_chain(object@chain[seq_len(i), j])
    }

    invisible(NULL)
}

##' @noRd
pmcmc_proposal <- function(object, i, theta_mean, covmat_emp,
                           scale_start = 100L, shape_start = 200L,
                           cooling = 0.999, max_scaling = 50) {
    n_accepted <- sum(object@chain[seq_len(i - 1), "accept"])
    n_pars <- length(object@pars)
    j <- seq(from = 5, by = 1, length.out = n_pars)
    theta <- object@chain[i - 1, j]

    if (runif(1) < object@adaptmix || i <= scale_start) {
        covmat <- diag((object@chain[1, j] / 10)^2 / n_pars, n_pars)
    } else if (n_accepted < shape_start) {
        scaling <- 1
        target <- ifelse(n_pars == 1L, 0.44, 0.234)

        for (k in seq(from = scale_start, to = i - 1L)) {
            l <- cooling^(k - scale_start)
            m <- mean(object@chain[seq_len(k), "accept"]) - target
            scaling <- min(scaling * exp(l * m), max_scaling)
        }

        covmat <- scaling^2 * diag((object@chain[1, j] / 10)^2 / n_pars, n_pars)
    } else {
        covmat <- 2.38^2 / n_pars * covmat_emp
    }

    theta_mean <- ((i - 1) * theta_mean + theta) / i
    covmat_emp <- ((i - 1) * covmat_emp + tcrossprod(theta - theta_mean)) / i

    ch <- chol(covmat, pivot = TRUE)
    Q <- ch[, order(attr(ch, "pivot"))]
    proposal <- as.numeric(theta + rnorm(n = n_pars) %*% Q)
    names(proposal) <- names(theta)

    list(theta = proposal,
         theta_mean = theta_mean,
         covmat_emp = covmat_emp)
}

covmat_empirical <- function(object, i) {
    n_pars <- length(object@pars)
    j <- seq(from = 5, by = 1, length.out = n_pars)
    covmat <- stats::cov(object@chain[seq_len(i), j, drop = FALSE])
    if (i == 1)
        covmat[, ] <- 0
    covmat
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
##' @template niter-param
##' @export
setMethod(
    "continue",
    signature(object = "SimInf_pmcmc"),
    function(object, niter, ...,
             verbose = getOption("verbose", FALSE)) {
        check_integer_arg(niter)
        niter <- as.integer(niter)
        if (any(length(niter) != 1L, any(niter <= 0L)))
            stop("'niter' must be an integer > 0.", call. = FALSE)

        iterations <- length(object) + seq_len(niter)
        object@chain <- setup_chain(object, niter)
        object@pf <- setup_pf(object, niter)

        if (iterations[1] == 1) {
            iterations <- iterations[-1]

            theta <- rpriors(object@priors)
            methods::slot(object@model, object@target) <-
                set_proposal(object, theta)

            object@pf[[1]] <- pfilter(object@model, object@obs_process,
                                      object@data, object@npart)

            pf <- object@pf[[1]]
            logLik <- pf@loglik
            logPrior <- dpriors(theta, object@priors)
            logPost <- logLik + logPrior
            accept <- 0
            object@chain[1, ] <- c(logPost, logLik, logPrior, accept, theta)
        }

        ## Continue from the last iteration in the chain.
        i <- iterations[1] - 1
        pf <- object@pf[[i]]
        logPost <- object@chain[i, "logPost"]
        logLik <- object@chain[i, "logLik"]
        logPrior <- object@chain[i, "logPrior"]
        j <- seq(from = 5, by = 1, length.out = length(object@pars))
        theta <- object@chain[i, j]
        theta_mean <- colMeans(object@chain[seq_len(i), j, drop = FALSE])
        covmat_emp <- covmat_empirical(object, i)

        for (i in iterations) {
            ## Proposal
            accept <- 0
            proposal <- pmcmc_proposal(object, i, theta_mean, covmat_emp)
            theta_mean <- proposal$theta_mean
            covmat_emp <- proposal$covmat_emp
            logPrior_prop <- dpriors(proposal$theta, object@priors)

            if (is.finite(logPrior_prop)) {
                methods::slot(object@model, object@target) <-
                    set_proposal(object, proposal$theta)

                pf_prop <- pfilter(object@model, object@obs_process,
                                   object@data, object@npart)
                logLik_prop <- pf_prop@loglik

                alpha <- exp(logLik_prop + logPrior_prop - logLik - logPrior)
                if (is.finite(alpha) && runif(1) < alpha) {
                    logLik <- logLik_prop
                    logPrior <- logPrior_prop
                    logPost <- logLik + logPrior
                    theta <- proposal$theta
                    pf <- pf_prop
                    accept <- 1
                }
            }

            ## Save current value of chain.
            object@chain[i, ] <- c(logPost, logLik, logPrior, accept, theta)
            object@pf[[i]] <- pf

            ## Report progress.
            pmcmc_progress(object, i, verbose)
        }

        object
    }
)

pmcmc_iterations <- function(x, start, end, thin) {
    check_integer_arg(start)
    start <- as.integer(start)
    if (any(length(start) != 1, any(start < 1)))
        stop("'start' must be an integer >= 1.", call. = FALSE)

    if (is.null(end))
        end <- length(x)
    check_integer_arg(end)
    end <- as.integer(end)
    if (any(length(end) != 1, any(end < start), any(end > length(x)))) {
        stop("'end' must be an integer between start and length(x).",
             call. = FALSE)
    }

    check_integer_arg(thin)
    thin <- as.integer(thin)
    if (any(length(thin) != 1, any(thin < 1)))
        stop("'thin' must be an integer >= 1.", call. = FALSE)

    seq(from = start, to = end, by = thin)
}
