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
##' @slot n_particles An integer with the number of particles (> 1) to
##'     use in the bootstrap particle filter.
##' @slot obs_process A \code{formula} or \code{function} determining
##'     the observation process.
##' @slot init_model An optional function that, if non-NULL, is
##'     applied before running each proposal. The function must accept
##'     one argument of type \code{SimInf_model} with the current
##'     model of the fitting process. This function can be useful to
##'     specify the initial state of \code{u0} or \code{v0} of the
##'     model before running a trajectory with proposed parameters.
##' @slot data A \code{data.frame} holding the time series data for
##'     the observation process.
##' @slot chain A matrix where each row contains \code{logPost},
##'     \code{logLik}, \code{logPrior}, \code{accept}, and the
##'     \code{parameters} for each iteration.
##' @slot pf List with the filtered trajectory from each iteration.
##' @slot covmat A named numeric \code{(npars x npars)} matrix with
##'     covariances to use as initial proposal matrix.
##' @slot adaptmix Mixing proportion for adaptive proposal.
##' @slot adaptive Controls when to start adaptive update.
##' @seealso \code{\link{pmcmc}} and \code{\link{continue_pmcmc}}.
##' @export
setClass(
    "SimInf_pmcmc",
    slots = c(model         = "SimInf_model",
              priors        = "data.frame",
              target        = "character",
              pars          = "integer",
              n_particles   = "integer",
              obs_process   = "ANY",
              init_model    = "ANY",
              data          = "data.frame",
              chain         = "matrix",
              pf            = "ANY",
              covmat        = "matrix",
              adaptmix      = "numeric",
              adaptive      = "integer")
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
        errors <- c(errors, "'adaptmix' must be a value >= 0 and <= 1.")
    }

    if (all(!identical(object@target, "gdata"),
            !identical(object@target, "ldata"))) {
        errors <- c(errors, "'target' must be 'gdata' or 'ldata'.")
    }

    if (all(!is.null(object@init_model),
            !is.function(object@init_model))) {
        errors <- c(errors, "'init_model' must be 'NULL' or a 'function'.")
    }

    if (length(errors))
        return(errors)
    TRUE
}

## Assign the validity method for the SimInf_pmcmc class.
setValidity("SimInf_pmcmc", valid_SimInf_pmcmc_object)

setAs(
    from = "SimInf_pmcmc",
    to = "data.frame",
    def = function(from) {
        ## Skip the first four columns in chain: 'logPost',
        ## 'logLik', 'logPrior', and 'accept'.
        j <- seq(from = 5, by = 1, length.out = n_pars(from))
        as.data.frame(from@chain[, j, drop = FALSE])
    }
)

##' Coerce to data frame
##'
##' @method as.data.frame SimInf_pmcmc
##'
##' @inheritParams base::as.data.frame
##' @export
as.data.frame.SimInf_pmcmc <- function(x, ...) {
    methods::as(x, "data.frame")
}

summary_chain <- function(chain) {
    qq <- do.call("rbind", apply(chain, 2, function(x) {
        cbind(t(quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)),
              Mean = mean(x, na.rm = TRUE),
              SD = sqrt(var(x, na.rm = TRUE)))
    }, simplify = FALSE))
    rownames(qq) <- colnames(chain)
    print.table(qq, digits = 3)
}

acceptance_ratio <- function(object) {
    mean(object@chain[, "accept"], na.rm = TRUE)
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
        cat(sprintf("Number of particles: %i\n", object@n_particles))
        cat(sprintf("Mixing proportion for adaptive proposal: %.2f\n",
                    object@adaptmix))

        if (length(object) > 0) {
            cat(sprintf("Acceptance ratio: %.3f\n", acceptance_ratio(object)))

            print_title(
                "Quantiles, mean and standard deviation for each variable")

            ## Skip first four columns in chain.
            j <- seq(from = 5, by = 1, length.out = n_pars(object))
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
        cat(sprintf("Number of particles: %i\n", object@n_particles))
        if (length(object) > 0)
            cat(sprintf("Acceptance ratio: %.3f\n", acceptance_ratio(object)))

        ## The model name
        cat(sprintf("Model: %s\n", as.character(class(object@model))))

        ## Nodes
        cat(sprintf("Number of nodes: %i\n", n_nodes(object)))

        summary_transitions(object@model)

        if (length(object) > 0) {
            print_title(
                "Quantiles, mean and standard deviation for each variable")

            ## Skip the first four columns in chain: 'logPost',
            ## 'logLik', 'logPrior', and 'accept'.
            j <- seq(from = 5, by = 1, length.out = n_pars(object))
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
##' @template n_particles-param
##' @template n_iterations-param
##' @param theta A named vector of initial values for the parameters
##'     of the model.  Default is \code{NULL}, and then these are
##'     sampled from the prior distribution(s).
##' @param covmat A named numeric \code{(npars x npars)} matrix with
##'     covariances to use as initial proposal matrix. If left
##'     unspecified then defaults to \code{diag((theta/10)^2/npars)}.
##' @param adaptmix Mixing proportion for adaptive proposal.  Must be
##'     a value between zero and one. Default is \code{adaptmix =
##'     0.05}.
##' @param adaptive Controls when to start adaptive update. Must be
##'     greater or equal to zero. If \code{adaptive=0}, then adaptive
##'     update is not performed. Default is \code{adaptive = 100}.
##' @template init_model-param
##' @param post_particle An optional function that, if non-NULL, is
##'     applied after each completed particle. The function must
##'     accept three arguments: 1) an object of \code{SimInf_pmcmc}
##'     with the current state of the fitting process, 2) an object
##'     \code{SimInf_pfilter} with the last particle and one filtered
##'     trajectory attached, and 3) an integer with the iteration in
##'     the fitting process. This function can be useful to, for
##'     example, monitor, save and inspect intermediate results.
##' @param record Record a filtered trajectory from the particle
##'     filter in each iteration. Default is \code{FALSE}.
##' @template verbose-param-pmcmc
##' @references
##'
##' \Andrieu2010
##'
##' \Roberts2009
##' @export
##' @seealso \code{\link{continue_pmcmc}}.
setGeneric(
    "pmcmc",
    signature = "model",
    function(model,
             obs_process,
             data,
             priors,
             n_particles,
             n_iterations,
             theta = NULL,
             covmat = NULL,
             adaptmix = 0.05,
             adaptive = 100,
             init_model = NULL,
             post_particle = NULL,
             record = FALSE,
             verbose = getOption("verbose", FALSE)) {
        standardGeneric("pmcmc")
    }
)

##' @rdname pmcmc
##' @export
setMethod(
    "pmcmc",
    signature(model = "SimInf_model"),
    function(model,
             obs_process,
             data,
             priors,
             n_particles,
             n_iterations,
             theta,
             covmat,
             adaptmix,
             adaptive,
             init_model,
             post_particle,
             record,
             verbose) {
        n_particles <- check_n_particles(n_particles)
        n_iterations <- check_n_iterations(n_iterations)
        adaptmix <- check_adaptmix(adaptmix)
        adaptive <- check_adaptive(adaptive)
        init_model <- check_init_model(init_model)
        post_particle <- check_post_particle(post_particle)
        pf <- check_record_pf(record)

        ## Match the 'priors' to parameters in 'ldata' or 'gdata'.
        priors <- parse_priors(priors)
        pars <- match_priors(model, priors)

        if (is.null(theta))
            theta <- rpriors(priors)
        if (!all(is.atomic(theta),
                 is.numeric(theta),
                 all(priors$parameter %in% names(theta)))) {
            stop("'theta' must be a vector with initial ",
                 "values for the parameters.",
                 call. = FALSE)
        }
        theta <- theta[priors$parameter]

        if (is.null(covmat)) {
            covmat <- diag(((theta / 10)^2) / length(theta),
                           nrow = length(theta))
            colnames(covmat) <- names(theta)
            rownames(covmat) <- names(theta)
        }

        object <- new("SimInf_pmcmc",
                      model = model,
                      priors = priors,
                      target = pars$target,
                      pars = pars$pars,
                      obs_process = obs_process,
                      init_model = init_model,
                      pf = pf,
                      data = data,
                      n_particles = n_particles,
                      covmat = covmat,
                      adaptmix = adaptmix,
                      adaptive = adaptive)

        object@chain <- setup_chain(object, 1L)
        object@pf <- setup_pf(object, 1L)

        methods::slot(object@model, object@target) <-
            set_proposal(object, theta)

        pf <- pfilter(object@model,
                      obs_process = object@obs_process,
                      data = object@data,
                      n_particles = object@n_particles,
                      init_model = object@init_model)

        logLik <- pf@loglik
        logPrior <- dpriors(theta, object@priors)
        logPost <- logLik + logPrior
        accept <- 0

        ## Save current value of chain.
        object@chain[1, ] <- c(logPost, logLik, logPrior, accept, theta)
        if (is.function(post_particle))
            post_particle(object, pf, 1)
        if (!is.null(object@pf))
            object@pf[[1]] <- pf

        n_iterations <- n_iterations - 1L
        if (n_iterations == 0)
            return(object)

        continue_pmcmc(object,
                       n_iterations = n_iterations,
                       post_particle = post_particle,
                       verbose = verbose)
    }
)

check_adaptive <- function(adaptive) {
    check_integer_arg(adaptive)
    adaptive <- as.integer(adaptive)
    if (any(length(adaptive) != 1L, any(adaptive < 0L)))
        stop("'adaptive' must be an integer >= 0.", call. = FALSE)
    adaptive
}

check_adaptmix <- function(adaptmix) {
    adaptmix <- as.numeric(adaptmix)
    if (any(length(adaptmix) != 1L,
            any(adaptmix <= 0),
            any(adaptmix >= 1))) {
        stop("'adaptmix' must be a value > 0 and < 1.", call. = FALSE)
    }
    adaptmix
}

check_init_model <- function(init_model) {
    if (!is.null(init_model))
        init_model <- match.fun(init_model)
    init_model
}

check_n_iterations <- function(n_iterations) {
    check_integer_arg(n_iterations)
    n_iterations <- as.integer(n_iterations)
    if (any(length(n_iterations) != 1L, any(n_iterations <= 0L)))
        stop("'n_iterations' must be an integer > 0.", call. = FALSE)
    n_iterations
}

check_post_particle <- function(post_particle) {
    if (!is.null(post_particle))
        post_particle <- match.fun(post_particle)
    post_particle
}

check_record_pf <- function(record) {
    if (isTRUE(record))
        return(list())
    NULL
}

is_empty_chain <- function(object) {
    isTRUE(length(object) == 0L)
}

setup_chain <- function(object, n_iterations) {
    m <- matrix(NA_real_,
                nrow = n_iterations,
                ncol = 4L + n_pars(object),
                dimnames = list(NULL, c("logPost", "logLik",
                                        "logPrior", "accept",
                                        object@priors$parameter)))

    if (is_empty_chain(object))
        return(m)
    rbind(object@chain, m)
}

setup_pf <- function(object, n_iterations) {
    if (is.null(object@pf))
        return(NULL)
    c(object@pf, vector("list", length = n_iterations))
}

set_proposal <- function(object, theta) {
    if (object@target == "gdata") {
        for (i in seq_len(n_pars(object))) {
            object@model@gdata[object@pars[i]] <- theta[i]
        }
    } else {
        for (i in seq_len(n_pars(object))) {
            object@model@ldata[object@pars[i], ] <- theta[i]
        }
    }

    methods::slot(object@model, object@target)
}

pmcmc_progress <- function(object, i, verbose) {
    if (!is.null(verbose) && isTRUE(i %% verbose == 0)) {
        print_title(sprintf(
            "Iteration: %i of %i. Time: %s. Acceptance ratio: %.3f",
            i, length(object), format(Sys.time(), "%T"),
            acceptance_ratio(object)))

        ## Skip columns logLik, logPrior and accept in the chain.
        j <- c(1, seq(from = 5, by = 1, length.out = n_pars(object)))
        summary_chain(object@chain[seq_len(i), j])
    }

    invisible(NULL)
}

n_pars <- function(x) {
    length(x@pars)
}

get_theta <- function(x, i) {
    j <- seq(from = 5, by = 1, length.out = n_pars(x))
    x@chain[i, j]
}

##' @noRd
pmcmc_proposal <- function(x, i, n_accepted, theta_mean, covmat_emp) {
    if (x@adaptive == 0L ||
        i <= x@adaptive ||
        n_accepted == 0 ||
        runif(1) < x@adaptmix) {
        covmat <- x@covmat
    } else {
        covmat <- 2.38^2 / n_pars(x) * covmat_emp
    }

    theta <- get_theta(x, i - 1)
    theta_mean <- ((i - 1) * theta_mean + theta) / i
    covmat_emp <- ((i - 1) * covmat_emp + tcrossprod(theta - theta_mean)) / i

    ch <- chol(covmat, pivot = TRUE)
    Q <- ch[, order(attr(ch, "pivot"))]
    proposal <- as.numeric(theta + rnorm(n = n_pars(x)) %*% Q)
    names(proposal) <- names(theta)

    list(theta = proposal,
         theta_mean = theta_mean,
         covmat_emp = covmat_emp)
}

covmat_empirical <- function(object, i) {
    j <- seq(from = 5, by = 1, length.out = n_pars(object))
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

get_verbose <- function(verbose) {
    if (isTRUE(verbose))
        return(100L)

    if (all(is.numeric(verbose),
            !anyNA(verbose),
            all(is_wholenumber(verbose)),
            all(verbose > 0),
            length(verbose) == 1L)) {
        return(as.integer(verbose))
    }

    NULL
}

##' Run more iterations of PMCMC
##'
##' @param object The \code{SimInf_pmcmc} object to continue from.
##' @template n_iterations-param
##' @param post_particle An optional function that, if non-NULL, is
##'     applied after each completed particle. The function must
##'     accept three arguments: 1) an object of \code{SimInf_pmcmc}
##'     with the current state of the fitting process, 2) an object
##'     \code{SimInf_pfilter} with the last particle and one filtered
##'     trajectory attached, and 3) an integer with the iteration in
##'     the fitting process. This function can be useful to, for
##'     example, monitor, save and inspect intermediate results.
##' @template verbose-param-pmcmc
##' @export
setGeneric(
    "continue_pmcmc",
    signature = "object",
    function(object,
             n_iterations,
             post_particle = NULL,
             verbose = getOption("verbose", FALSE)) {
        standardGeneric("continue_pmcmc")
    }
)

##' @rdname continue_pmcmc
##' @export
setMethod(
    "continue_pmcmc",
    signature(object = "SimInf_pmcmc"),
    function(object,
             n_iterations,
             post_particle,
             verbose) {
        methods::validObject(object)

        n_iterations <- check_n_iterations(n_iterations)
        post_particle <- check_post_particle(post_particle)
        verbose <- get_verbose(verbose)
        iterations <- length(object) + seq_len(n_iterations)
        object@chain <- setup_chain(object, n_iterations)
        object@pf <- setup_pf(object, n_iterations)

        ## Continue from the last iteration in the chain.
        i <- iterations[1] - 1
        if (!is.null(object@pf))
            pf <- object@pf[[i]]
        n_accepted <- sum(object@chain[seq_len(i), "accept"])
        logPost <- object@chain[i, "logPost"]
        logLik <- object@chain[i, "logLik"]
        logPrior <- object@chain[i, "logPrior"]
        j <- seq(from = 5, by = 1, length.out = n_pars(object))
        theta <- object@chain[i, j]
        theta_mean <- colMeans(object@chain[seq_len(i), j, drop = FALSE])
        covmat_emp <- covmat_empirical(object, i)

        for (i in iterations) {
            ## Proposal
            accept <- 0
            proposal <- pmcmc_proposal(x = object,
                                       i = i,
                                       n_accepted = n_accepted,
                                       theta_mean = theta_mean,
                                       covmat_emp = covmat_emp)
            theta_mean <- proposal$theta_mean
            covmat_emp <- proposal$covmat_emp
            logPrior_prop <- dpriors(proposal$theta, object@priors)

            if (is.finite(logPrior_prop)) {
                methods::slot(object@model, object@target) <-
                    set_proposal(object, proposal$theta)

                pf_prop <- pfilter(model = object@model,
                                   obs_process = object@obs_process,
                                   data = object@data,
                                   n_particles = object@n_particles,
                                   init_model = object@init_model)
                logLik_prop <- pf_prop@loglik

                alpha <- exp(logLik_prop + logPrior_prop - logLik - logPrior)
                if (is.finite(alpha) && runif(1) < alpha) {
                    logLik <- logLik_prop
                    logPrior <- logPrior_prop
                    logPost <- logLik + logPrior
                    theta <- proposal$theta
                    pf <- pf_prop
                    accept <- 1
                    n_accepted <- n_accepted + 1
                }
            }

            ## Save current value of chain.
            object@chain[i, ] <- c(logPost, logLik, logPrior, accept, theta)
            if (is.function(post_particle))
                post_particle(object, pf, i)
            if (!is.null(object@pf))
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
