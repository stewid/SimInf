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

##' Class \code{SimInf_if2}
##'
##' Storage class for the results of iterated filtering (IF2) on a
##' \code{SimInf_model}, see \code{\link{if2}}.
##'
##' @slot model The \code{SimInf_model} object to estimate parameters
##'     in.
##' @slot n_particles An integer with the number of particles in the
##'     parameter swarm.
##' @slot cooling_fraction A numeric scalar with the fraction of the
##'     perturbation scale that remains after 50 iterations. The
##'     perturbation scales are multiplied by
##'     \code{cooling_fraction^(i/50)} at iteration \code{i}.
##' @slot loglik A numeric vector with the estimated log-likelihood of
##'     the perturbed model for each iteration of the iterated
##'     filtering. The final estimate is available with
##'     \code{\link[=logLik,SimInf_if2-method]{logLik}}. Note that the
##'     estimates are inflated close to the maximum and are not
##'     suitable for likelihood-based inference on the focal model.
##' @slot ess A matrix with the effective sample size, with one row
##'     per observation time and one column per iteration. Entries are
##'     \code{NA} for intervals that were not reached in an iteration,
##'     e.g. if the swarm degenerated.
##' @include SimInf_model.R
##' @export
setClass(
    "SimInf_if2",
    slots = c(model = "SimInf_model",
              n_particles = "integer",
              cooling_fraction = "numeric",
              loglik = "numeric",
              ess = "matrix")
)

##' Brief summary of a \code{SimInf_if2} object
##'
##' @param object The \code{SimInf_if2} object.
##' @return \code{invisible(object)}.
##' @export
setMethod(
    "show",
    signature(object = "SimInf_if2"),
    function(object) {
        cat(sprintf("Number of particles: %i\n", object@n_particles))
        cat(sprintf("Log-likelihood: %f\n", logLik(object)))

        invisible(object)
    }
)

##' Detailed summary of a \code{SimInf_if2} object
##'
##' @param object The \code{SimInf_if2} object
##' @param ... Not used.
##' @return None (invisible 'NULL').
##' @export
setMethod(
    "summary",
    signature(object = "SimInf_if2"),
    function(object, ...) {
        cat("Iterated filtering\n")
        cat("------------------\n")
        cat(sprintf("Number of particles: %i\n", object@n_particles))
        cat(sprintf("Cooling fraction: %.2f\n", object@cooling_fraction))
        cat(sprintf("Log-likelihood: %f\n", logLik(object)))

        invisible(NULL)
    }
)

##' Run iterated filtering on a model.
##' @noRd
if2_internal <- function(model,
                         events,
                         obs_process,
                         data,
                         n_particles,
                         tspan,
                         cooling,
                         init_model,
                         post_particle) {
    Nd <- nrow(model@v0)
    Ntspan <- nrow(tspan)
    ess <- rep(NA_real_, Ntspan)
    loglik <- 0
    w <- numeric(n_particles)

    ## Keep track of the states in the compartments and the continuous
    ## states in every particle. Only the states from the latest
    ## observation time are needed to rebuild u0 and v0 for the coming
    ## interval.
    if (is.function(init_model)) {
        U <- integer(n_particles * n_nodes(model) * n_compartments(model))
        V <- numeric(n_particles * n_nodes(model) * Nd)

        ## Loop over the particles and initialise each model.
        for (p in seq_len(n_particles)) {
            ## Initialise the model.
            m <- init_model(model)

            ## Guard against 'init_model' changing the structure of
            ## the model: the particle blocks in 'U' and 'V' are
            ## filled with the template layout, so 'u0' and 'v0' must
            ## keep the dimensions of the template.
            if (!identical(dim(m@u0), dim(model@u0)) ||
                !identical(dim(m@v0), dim(model@v0))) {
                stop("'init_model' must not change the dimensions ",
                     "of 'u0' or 'v0'.", call. = FALSE)
            }

            ## Save the initial u0 state.
            u_i <- seq.int(
                from = (p - 1L) * n_nodes(model) * n_compartments(model) + 1L,
                length.out = n_nodes(model) * n_compartments(model))
            U[u_i] <- as.integer(m@u0)

            ## Save the initial v0 state.
            v_i <- seq.int(
                from = (p - 1L) * n_nodes(model) * Nd + 1L,
                length.out = n_nodes(model) * Nd)
            V[v_i] <- as.numeric(m@v0)
        }
    } else {
        U <- rep(as.integer(model@u0), n_particles)
        V <- rep(as.numeric(model@v0), n_particles)
    }

    ## Loop over time series.
    a <- seq_len(n_particles)
    m <- model
    m@replicates <- n_particles
    for (i in seq_len(Ntspan)) {
        if (is.na(tspan[i, 1L])) {
            m@tspan <- tspan[i, 2L]
        } else {
            m@tspan <- tspan[i, 1:2]
        }

        ## Initialise events for the interval.
        if (!is.null(events))
            m@events <- events[[i]]

        ## Initialise the model.
        u_i <- seq_len(n_nodes(model) * n_compartments(model)) +
            rep((a - 1L) * n_nodes(model) * n_compartments(model),
                each = n_nodes(model) * n_compartments(model))
        m@u0 <- matrix(
            data = U[u_i],
            nrow = n_compartments(model),
            dimnames = dimnames(m@u0))

        v_i <- seq_len(n_nodes(model) * Nd) +
            rep((a - 1L) * n_nodes(model) * Nd,
                each = n_nodes(model) * Nd)
        m@v0 <- matrix(
            data = V[v_i],
            nrow = nrow(m@v0),
            dimnames = dimnames(m@v0))

        ## Propagate the model.
        x <- run(m)
        if (length(x@tspan) > 1L) {
            ## Extract the final state in each replicate. The first
            ## replicate is in column 1 and 2, the second replicate in
            ## column 3 and 4, et cetera.
            x@tspan <- x@tspan[2L]
            j <- seq.int(from = 2L, by = 2L, length.out = n_particles)
            x@U <- x@U[, j, drop = FALSE]
            x@V <- x@V[, j, drop = FALSE]
        }

        ## Save the final state of all particles. Writing into the
        ## preallocated buffers keeps them as flat vectors, with the
        ## values placed in column-major order, which matches the
        ## particle-block layout used when rebuilding u0 and v0.
        U[] <- x@U
        V[] <- x@V

        ## Set the weight for the particle.
        if (is.function(obs_process)) {
            w <- obs_process(x, data[[i]])
        } else {
            e <- new.env(parent = baseenv())

            assign(x = obs_process$par,
                   value = data[[i]][, obs_process$par_i],
                   pos = e)

            for (j in seq_along(obs_process$slots)) {
                s <- obs_process$slots[[j]]$slot
                v <- slot(x, s)[obs_process$slots[[j]]$i, ]
                n <- obs_process$slots[[j]]$name
                assign(x = n, value = v, pos = e)
            }

            expr <- obs_process$expr
            e$expr <- expr
            w <- evalq(eval(parse(text = expr)), envir = e)
        }

        ## The observation process must return a numeric vector with
        ## one entry per particle. -Inf is a valid outcome: the
        ## observation is impossible under that particle, and the
        ## weight normalisation assigns it zero weight. +Inf/NA/NaN
        ## are invalid.
        if (!identical(length(w), n_particles) || !is.vector(w, "numeric")) {
            stop("Invalid observation process: weights must be ",
                 "a numeric vector of length 'n_particles'.",
                 call. = FALSE)
        }
        if (anyNA(w)) {
            stop("Invalid observation process: weights must not be NA or NaN.",
                 call. = FALSE)
        }
        if (any(w == Inf)) {
            stop("Invalid observation process: weights must not be +Inf.",
                 call. = FALSE)
        }

        max_w <- max(w)
        if (!is.finite(max_w)) {
            ## Signal a degenerate swarm to the outer loop.
            return(list(ok = FALSE,
                        loglik = NA_real_,
                        ess = ess))
        }

        w <- exp(w - max_w)
        sum_w <- sum(w)
        loglik <- loglik + max_w + log(sum_w) - log(n_particles)
        w <- w / sum_w
        ess[i] <- 1 / sum(w^2)

        ## Systematic resampling. Sample the ancestor indices for the
        ## next generation of particles from the normalised weights.
        if (i < Ntspan) {
            a[] <- .Call(SimInf_systematic_resampling, w)
        }
    }

    list(ok = TRUE,
         loglik = loglik,
         ess = ess)
}

##' Run iterated filtering on a model.
##' @noRd
if2_ldata <- function(model,
                      events,
                      obs_process,
                      data,
                      tspan,
                      cooling,
                      init_model,
                      post_particle) {
    Nd <- nrow(model@v0)
    Ntspan <- nrow(tspan)
    ess <- rep(NA_real_, Ntspan)
    loglik <- 0
    w <- numeric(n_nodes(model))

    ## Optionally initialise the initial state of every particle with
    ## the 'init_model' callback. The callback is called once per
    ## iteration with the expanded model, and only 'u0' and 'v0' are
    ## extracted from the returned model: writes to any other slot,
    ## including 'ldata', 'gdata' and 'events', are discarded. Only
    ## the states from the latest observation time are needed to
    ## rebuild u0 and v0 for the coming interval.
    if (is.function(init_model)) {
        ## Initialise the model.
        m <- init_model(model)
        if (!methods::is(m, "SimInf_model")) {
            stop("'init_model' must return a 'SimInf_model' object.",
                 call. = FALSE)
        }

        ## Guard against 'init_model' changing the structure of the
        ## model: each node in the expanded model is one particle, so
        ## 'u0' and 'v0' must keep the dimensions of the template for
        ## the per-particle states to stay aligned with the
        ## per-particle parameters in 'ldata'.
        if (!identical(dim(m@u0), dim(model@u0)) ||
            !identical(dim(m@v0), dim(model@v0))) {
            stop("'init_model' must not change the dimensions ",
                 "of 'u0' or 'v0'.", call. = FALSE)
        }

        model@u0 <- matrix(
            data = as.integer(m@u0),
            nrow = n_compartments(model),
            dimnames = dimnames(model@u0))

        model@v0 <- matrix(
            data = as.numeric(m@v0),
            nrow = nrow(model@v0),
            dimnames = dimnames(model@v0))
    }

    ## Loop over time series.
    a <- seq_len(n_particles)
    m <- model
    m@replicates <- n_particles
    for (i in seq_len(Ntspan)) {
        if (is.na(tspan[i, 1L])) {
            m@tspan <- tspan[i, 2L]
        } else {
            m@tspan <- tspan[i, 1:2]
        }

        ## Initialise events for the interval.
        if (!is.null(events))
            m@events <- events[[i]]

        ## Initialise the model.
        u_i <- seq_len(n_nodes(model) * n_compartments(model)) +
            rep((a - 1L) * n_nodes(model) * n_compartments(model),
                each = n_nodes(model) * n_compartments(model))
        m@u0 <- matrix(
            data = U[u_i],
            nrow = n_compartments(model),
            dimnames = dimnames(m@u0))

        v_i <- seq_len(n_nodes(model) * Nd) +
            rep((a - 1L) * n_nodes(model) * Nd,
                each = n_nodes(model) * Nd)
        m@v0 <- matrix(
            data = V[v_i],
            nrow = nrow(m@v0),
            dimnames = dimnames(m@v0))

        ## Propagate the model.
        x <- run(m)
        if (length(x@tspan) > 1L) {
            ## Extract the final state in each replicate. The first
            ## replicate is in column 1 and 2, the second replicate in
            ## column 3 and 4, et cetera.
            x@tspan <- x@tspan[2L]
            j <- seq.int(from = 2L, by = 2L, length.out = n_particles)
            x@U <- x@U[, j, drop = FALSE]
            x@V <- x@V[, j, drop = FALSE]
        }

        ## Save the final state of all particles. Writing into the
        ## preallocated buffers keeps them as flat vectors, with the
        ## values placed in column-major order, which matches the
        ## particle-block layout used when rebuilding u0 and v0.
        U[] <- x@U
        V[] <- x@V

        ## Set the weight for the particle.
        if (is.function(obs_process)) {
            w <- obs_process(x, data[[i]])
        } else {
            e <- new.env(parent = baseenv())

            assign(x = obs_process$par,
                   value = data[[i]][, obs_process$par_i],
                   pos = e)

            for (j in seq_along(obs_process$slots)) {
                s <- obs_process$slots[[j]]$slot
                v <- slot(x, s)[obs_process$slots[[j]]$i, ]
                n <- obs_process$slots[[j]]$name
                assign(x = n, value = v, pos = e)
            }

            expr <- obs_process$expr
            e$expr <- expr
            w <- evalq(eval(parse(text = expr)), envir = e)
        }

        ## The observation process must return a numeric vector with
        ## one entry per particle. -Inf is a valid outcome: the
        ## observation is impossible under that particle, and the
        ## weight normalisation assigns it zero weight. +Inf/NA/NaN
        ## are invalid.
        if (!identical(length(w), n_particles) || !is.vector(w, "numeric")) {
            stop("Invalid observation process: weights must be ",
                 "a numeric vector of length 'n_particles'.",
                 call. = FALSE)
        }
        if (anyNA(w)) {
            stop("Invalid observation process: weights must not be NA or NaN.",
                 call. = FALSE)
        }
        if (any(w == Inf)) {
            stop("Invalid observation process: weights must not be +Inf.",
                 call. = FALSE)
        }

        max_w <- max(w)
        if (!is.finite(max_w)) {
            ## Signal a degenerate swarm to the outer loop.
            return(list(ok = FALSE,
                        loglik = NA_real_,
                        ess = ess))
        }

        w <- exp(w - max_w)
        sum_w <- sum(w)
        loglik <- loglik + max_w + log(sum_w) - log(n_particles)
        w <- w / sum_w
        ess[i] <- 1 / sum(w^2)

        ## Systematic resampling. Sample the ancestor indices for the
        ## next generation of particles from the normalised weights.
        if (i < Ntspan) {
            a[] <- .Call(SimInf_systematic_resampling, w)
        }
    }

    list(ok = TRUE,
         loglik = loglik,
         ess = ess)
}

##' Iterated filtering algorithm
##'
##' @param model The model to simulate data from.
##' @template obs_process-param
##' @template data-param
##' @template n_particles-param
##' @param n_iterations An integer specifying the number of filtering
##'     iterations to perform.
##' @template verbose-param
##' @references
##'
##' \Ionides2015
##' @export
## nolint start: brace_linter
setGeneric(
    "if2",
    signature = "model",
    function(model,
             obs_process,
             data,
             n_particles,
             n_iterations,
             cooling_fraction = 0.5,
             init_model = NULL,
             post_particle = NULL,
             verbose = FALSE)
        standardGeneric("if2")
)
## nolint end

##' @rdname if2
##' @export
setMethod(
    "if2",
    signature(model = "SimInf_model"),
    function(model,
             obs_process,
             data,
             n_particles,
             n_iterations,
             cooling_fraction,
             init_model,
             post_particle,
             verbose) {
        if (isFALSE(identical(dim(model@U_sparse), c(0L, 0L))) ||
            isFALSE(identical(dim(model@V_sparse), c(0L, 0L)))) {
            stop("'if2' cannot run a model with a sparse result matrix.",
                 call. = FALSE)
        }

        n_particles <- check_n_particles(n_particles)
        n_iterations <- check_n_iterations(n_iterations, FALSE)
        cooling_fraction <- check_cooling_fraction(cooling_fraction)
        cooling_rate <- cooling_fraction^(1 / 50)
        data <- split_data(model, data)
        tspan <- split_tspan(model, data)
        model@tspan <- tspan[, 2]
        events <- split_events(model@events, tspan[, 2] - 1)
        verbose <- check_verbose(verbose)
        obs_process <- pfilter_obs_process(model, obs_process, data)
        init_model <- check_init_model(init_model)
        post_particle <- check_post_particle(post_particle)

        loglik <- rep(NA_real_, n_iterations)
        ess <- matrix(NA_real_, nrow = nrow(tspan), ncol = n_iterations)
        for (i in seq_len(n_iterations)) {
            cooling <- cooling_rate^(i - 1L)

            out <- if2_internal(model = model,
                                events = events,
                                obs_process = obs_process,
                                data = data,
                                n_particles = n_particles,
                                tspan = tspan,
                                cooling = cooling,
                                init_model = init_model,
                                post_particle = post_particle)

            loglik[i] <- out$loglik
            ess[, i] <- out$ess
            if (!isTRUE(out$ok)) {
                warning(sprintf(
                    "Degenerate swarm at iteration %i at interval %i",
                    i, min(which(is.na(out$ess)))),
                    call. = FALSE)
            }

            if (isTRUE(i %% verbose == 0L)) {
                message(sprintf(
                    "Iteration: %i of %i. Time: %s. Loglik: %.3f",
                    i, n_iterations, format(Sys.time(), "%T"),
                    loglik[i]))
            }
        }

        methods::new("SimInf_if2",
                     model = model,
                     n_particles = n_particles,
                     cooling_fraction = cooling_fraction,
                     loglik = loglik,
                     ess = ess)
    }
)

##' Extract the log-likelihood
##'
##' Extract the estimated log-likelihood from a \code{SimInf_if2}
##' object.
##'
##' @details The log-likelihood is recorded for each iteration of the
##'     iterated filtering, and this method returns the estimate from
##'     the final iteration. Note that iterated filtering estimates
##'     the likelihood of the perturbed model, which is typically
##'     inflated close to the maximum, so the returned value is not
##'     suitable for likelihood-based inference on the focal model. To
##'     obtain a log-likelihood estimate at the parameter values in
##'     the result, run
##'     \code{\link[=pfilter,SimInf_model-method]{pfilter}} on the
##'     model with the estimated parameters.
##'
##' @param object The \code{\linkS4class{SimInf_if2}} object.
##' @return The estimated log-likelihood, returned as a numeric
##'     scalar. Or \code{NA} if the final iteration was degenerate.
##' @seealso \code{\link{if2}} for running iterated filtering and
##'     creating \code{SimInf_if2} objects,
##'     \code{\linkS4class{SimInf_if2}} for the class definition.
##' @export
setMethod(
    "logLik",
    signature(object = "SimInf_if2"),
    function(object) {
        tail(object@loglik, n = 1)
    }
)
