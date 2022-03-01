## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2022 Stefan Widgren
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

##' Class \code{"SimInf_abc"}
##'
##' @slot model The \code{SimInf_model} object to estimate parameters
##'     in.
##' @template priors-slot
##' @slot target Character vector (\code{gdata} or \code{ldata}) that
##'     determines if the ABC-SMC method estimates parameters in
##'     \code{model@@gdata} or in \code{model@@ldata}.
##' @slot pars Index to the parameters in \code{target}.
##' @slot npart The number of particles in each generation.
##' @slot nprop An integer vector with the number of simulated
##'     proposals in each generation.
##' @slot fn A function for calculating the summary statistics for the
##'     simulated trajectory and determine for each particle if it
##'     should be accepted (\code{TRUE}) or rejected (\code{FALSE}).
##'     The first argument in \code{fn} is the simulated model
##'     containing one trajectory.  The second argument to \code{fn}
##'     is an integer with the \code{generation} of the particles.
##'     The function should return a logical vector with one value for
##'     each particle in the simulated model.
##' @slot tolerance A numeric matrix (number of summary statistics X
##'     number of generations) where each column contains the
##'     tolerances for a generation and each row contains a sequence
##'     of gradually decreasing tolerances.
##' @slot x A list where each item is a \code{matrix} with the
##'     accepted particles in each generation. Each column is one
##'     particle.
##' @slot w A list where each item is a vector with the weights for
##'     the particles \code{x} in the corresponding generation.
##' @slot distance A list where each item is a numeric matrix (number
##'     of summary statistics X number of particles) with the distance
##'     for the particles \code{x} in the corresponding generation.
##'     Each column contains the distance for a particle and each row
##'     contains the distance for a summary statistics.
##' @slot ess A numeric vector with the effective sample size (ESS) in
##'     each generation. Effective sample size is computed as
##'     \deqn{\left(\sum_{i=1}^N\!(w_{g}^{(i)})^2\right)^{-1},}{1 /
##'     (sum(w_ig^2)),} where \eqn{w_{g}^{(i)}}{w_ig} is the
##'     normalized weight of particle \eqn{i} in generation \eqn{g}.
##' @seealso \code{\link{abc}} and \code{\link{continue}}.
##' @export
setClass(
    "SimInf_abc",
    slots = c(model     = "SimInf_model",
              priors    = "data.frame",
              target    = "character",
              pars      = "integer",
              npart     = "integer",
              nprop     = "integer",
              fn        = "function",
              tolerance = "matrix",
              x         = "list",
              w         = "list",
              distance  = "list",
              ess       = "numeric")
)

setAs(
    from = "SimInf_abc",
    to = "data.frame",
    def = function(from) {
        do.call("rbind", lapply(seq_len(length(from@x)), function(i) {
            cbind(generation = i,
                  weight = from@w[[i]],
                  as.data.frame(t(from@x[[i]])))
        }))
    }
)

##' Coerce to data frame
##'
##' @method as.data.frame SimInf_abc
##'
##' @inheritParams base::as.data.frame
##' @export
as.data.frame.SimInf_abc <- function(x, ...) {
    as(x, "data.frame")
}

summary_particles <- function(object, i) {
    str <- sprintf("Generation %i:", i)
    cat(sprintf("\n%s\n", str))
    cat(sprintf("%s\n", paste0(rep("-", nchar(str)), collapse = "")))
    cat(sprintf(" Accrate: %.2e\n", object@npart / object@nprop[i]))
    cat(sprintf(" ESS: %.2e\n\n", object@ess[i]))
    summary_matrix(object@x[[i]])
}

##' Print summary of a \code{SimInf_abc} object
##'
##' @param object The \code{SimInf_abc} object.
##' @return \code{invisible(object)}.
##' @export
##' @importFrom methods show
setMethod(
    "show",
    signature(object = "SimInf_abc"),
    function(object) {
        cat(sprintf("Number of particles per generation: %i\n",
                    object@npart))
        cat(sprintf("Number of generations: %i\n",
                    length(object@x)))

        if (length(object@x))
            summary_particles(object, length(object@x))

        invisible(object)
    }
)

##' Detailed summary of a \code{SimInf_abc} object
##'
##' @param object The \code{SimInf_abc} object
##' @param ... Additional arguments affecting the summary produced.
##' @return None (invisible 'NULL').
##' @include SimInf_model.R
##' @export
setMethod(
    "summary",
    signature(object = "SimInf_abc"),
    function(object, ...) {
        cat(sprintf("Number of particles per generation: %i\n",
                    object@npart))
        cat(sprintf("Number of generations: %i\n",
                    length(object@x)))

        for (i in seq_len(length(object@x))) {
            summary_particles(object, i)
        }

        invisible(NULL)
    }
)

##' Generate replicates of the first node in the model.
##'
##' Replicate the node specific matrices 'u0', 'v0' and 'ldata' in the
##' first node. Additionally, replicate any events.
##' @param model the model to replicate.
##' @param n the number of replicates.
##' @param n_events the number of events for the first node in the
##'     model.
##' @return A modified model object
##' @noRd
replicate_first_node <- function(model, n, n_events) {
    if (dim(model@u0)[1] > 0)
        model@u0 <- model@u0[, rep(1, n), drop = FALSE]
    if (dim(model@v0)[1] > 0)
        model@v0 <- model@v0[, rep(1, n), drop = FALSE]
    if (dim(model@ldata)[1] > 0)
        model@ldata <- model@ldata[, rep(1, n), drop = FALSE]

    if (n_events > 0) {
        ## Replicate the events in the first node and add an offset to
        ## the node vector. The offset is not added to 'dest' since
        ## there are no external transfer events.
        i <- seq_len(n_events)
        offset <- rep(seq_len(n) - 1L, each = n_events)
        model@events@event <- rep(model@events@event[i], n)
        model@events@time <- rep(model@events@time[i], n)
        model@events@node <- rep(model@events@node[i], n) + offset
        model@events@dest <- rep(model@events@dest[i], n)
        model@events@n <- rep(model@events@n[i], n)
        model@events@proportion <- rep(model@events@proportion[i], n)
        model@events@select <- rep(model@events@select[i], n)
        model@events@shift <- rep(model@events@shift[i], n)
    }

    model
}

proposal_covariance <- function(x) {
    if (is.null(x))
        return(NULL)
    cov(t(x)) * 2
}

n_particles <- function(x) {
    if (is.null(x))
        return(0)
    ncol(x)
}

abc_tolerance <- function(tolerance, epsilon, ngen) {
    if (!is.numeric(tolerance))
        stop("'tolerance' must have non-negative values.", call. = FALSE)

    if (!is.matrix(tolerance))
        dim(tolerance) <- c(1L, length(tolerance))

    if (is.integer(tolerance))
        storage.mode(tolerance) <- "double"

    if (!identical(ncol(tolerance), ngen))
        stop("'tolerance' must have 'ngen' columns.", call. = FALSE)

    if (any(is.na(tolerance)) || any(tolerance < 0))
        stop("'tolerance' must have non-negative values.", call. = FALSE)

    tol <- NULL
    if (nrow(epsilon) > 0) {
        ## Check that the number of summary statistics is the same as
        ## in previous generations.
        if (nrow(tolerance) != nrow(epsilon))
            stop("Invalid dimension of 'tolerance'.", call. = FALSE)
        tol <- epsilon
    }

    ## Check that tolerance is a decreasing vector for each summary
    ## statistics.
    if (any(apply(cbind(tol, tolerance), 1, diff) >= 0))
        stop("'tolerance' must be a decreasing vector.", call. = FALSE)

    tolerance
}

abc_progress <- function(t0, t1, x, w, npart, nprop) {
    cat(sprintf("\n\n  accrate = %.2e, ESS = %.2e time = %.2f secs\n\n",
                npart / nprop, 1 / sum(w^2), (t1 - t0)[3]))
    summary_matrix(x)
}

##' Check the result from the ABC distance function.
##' @noRd
abc_distance <- function(distance, n) {
    if (!is.numeric(distance)) {
        stop("The result from the ABC distance function must be numeric.",
             call. = FALSE)
    }

    if (!is.matrix(distance))
        dim(distance) <- c(1L, length(distance))

    if (is.integer(distance))
        storage.mode(distance) <- "double"

    if (!identical(ncol(distance), n)) {
        stop("Invalid dimension of the result from the ABC distance function.",
             call. = FALSE)
    }

    if (any(is.na(distance)) || any(distance < 0)) {
        stop("The result from the ABC distance function must be non-negative.",
             call. = FALSE)
    }

    distance
}

##' Determine which particles to accept
##'
##' @param distance a numeric matrix (number of summary statistics X
##'     number of particles) with the distance for the particles. Each
##'     column contains the distance for a particle and each row
##'     contains the distance for a summary statistics.
##' @param tolerance a numeric vector with the tolerance for each
##'     summary statistics.
##' @return a logical vector of length ncol(distance) with TRUE for
##'     the particles to accept, else FALSE.
##' @noRd
abc_accept <- function(distance, tolerance) {
    if (!identical(nrow(distance), length(tolerance))) {
        stop("Mismatch between the number of summary statistics and tolerance.",
             call. = FALSE)
    }

    colSums(distance <= tolerance) == length(tolerance)
}

##' @importFrom utils setTxtProgressBar
##' @importFrom utils txtProgressBar
##' @noRd
abc_gdata <- function(model, pars, priors, npart, fn, generation,
                      tolerance, x, w, verbose, ...) {
    if (isTRUE(verbose)) {
        cat("\nGeneration", generation, "...\n")
        pb <- txtProgressBar(min = 0, max = npart, style = 3)
        t0 <- proc.time()
    }

    xx <- NULL
    ancestor <- NULL
    nprop <- 0L
    sigma <- proposal_covariance(x)

    while (n_particles(xx) < npart) {
        proposals <- .Call(SimInf_abc_proposals, priors$parameter,
                           priors$distribution, priors$p1, priors$p2,
                           1L, x, w, sigma)
        for (i in seq_len(nrow(proposals))) {
            model@gdata[pars[i]] <- proposals[i, 1]
        }

        d <- abc_distance(fn(run(model), generation = generation, ...), 1L)
        accept <- abc_accept(d, tolerance)
        nprop <- nprop + 1L
        if (isTRUE(accept)) {
            ## Collect accepted particle
            xx <- cbind(xx, as.matrix(model@gdata)[pars, 1, drop = FALSE])
            ancestor <- c(ancestor, attr(proposals, "ancestor")[1])
        }

        ## Report progress.
        if (isTRUE(verbose))
            setTxtProgressBar(pb, n_particles(xx))
    }

    ## Calculate weights.
    ww <- .Call(SimInf_abc_weights, priors$distribution, priors$p1,
                priors$p2, x[, ancestor], xx, w, sigma)

    ## Report progress.
    if (isTRUE(verbose))
        abc_progress(t0, proc.time(), xx, ww, npart, nprop)

    list(x = xx, w = ww, nprop = nprop)
}

##' @importFrom utils setTxtProgressBar
##' @importFrom utils txtProgressBar
##' @noRd
abc_ldata <- function(model, pars, priors, npart, fn, generation,
                      tolerance, x, w, verbose, ...) {
    ## Let each node represents one particle. Replicate the first node
    ## to run multiple particles simultaneously. Start with 10 x
    ## 'npart' and then increase the number adaptively based on the
    ## acceptance rate.
    n <- as.integer(10 * npart)
    n_events <- length(model@events@event)
    model <- replicate_first_node(model, n, n_events)

    if (isTRUE(verbose)) {
        cat("\nGeneration", generation, "...\n")
        pb <- txtProgressBar(min = 0, max = npart, style = 3)
        t0 <- proc.time()
    }

    xx <- NULL
    ancestor <- NULL
    nprop <- 0L
    sigma <- proposal_covariance(x)

    while (n_particles(xx) < npart) {
        if (all(n < 1e5L, nprop > 2L * n)) {
            ## Increase the number of particles that are simulated in
            ## each trajectory.
            n <- min(1e5L, n * 2L)
            model <- replicate_first_node(model, n, n_events)
        }

        proposals <- .Call(SimInf_abc_proposals, priors$parameter,
                           priors$distribution, priors$p1, priors$p2,
                           n, x, w, sigma)
        for (i in seq_len(nrow(proposals))) {
            model@ldata[pars[i], ] <- proposals[i, ]
        }

        d <- abc_distance(fn(run(model), generation = generation, ...), n)
        accept <- abc_accept(d, tolerance)

        ## Collect accepted particles making sure not to collect more
        ## than 'npart'.
        i <- cumsum(accept) + n_particles(xx)
        i <- which(i == npart)
        j <- which(accept)
        if (length(i)) {
            i <- min(i)
            j <- j[j <= i]
            nprop <- nprop + i
            xx <- cbind(xx, model@ldata[pars, j, drop = FALSE])
            ancestor <- c(ancestor, attr(proposals, "ancestor")[j])
        } else {
            nprop <- nprop + n
            if (length(j)) {
                xx <- cbind(xx, model@ldata[pars, j, drop = FALSE])
                ancestor <- c(ancestor, attr(proposals, "ancestor")[j])
            }
        }

        ## Report progress.
        if (isTRUE(verbose))
            setTxtProgressBar(pb, n_particles(xx))
    }

    ## Calculate weights.
    ww <- .Call(SimInf_abc_weights, priors$distribution, priors$p1,
                priors$p2, x[, ancestor], xx, w, sigma)

    ## Report progress.
    if (isTRUE(verbose))
        abc_progress(t0, proc.time(), xx, ww, npart, nprop)

    list(x = xx, w = ww, nprop = nprop)
}

##' Approximate Bayesian computation
##'
##' @param model The model to generate data from.
##' @template priors-param
##' @param ngen The number of generations of ABC-SMC to run.
##' @param npart An integer specifying the number of particles.
##' @param fn A function for calculating the summary statistics for a
##'     simulated trajectory. For each particle, the function must
##'     determine the distance and return that information. The first
##'     argument passed to the \code{fn} function is the result from a
##'     \code{run} of the model and it contains one trajectory. The
##'     second argument to \code{fn} is an integer with the
##'     \code{generation} of the particle(s). Depending on the
##'     underlying model structure, data for one or more particles
##'     have been generated in each call to \code{fn}. If the
##'     \code{model} only contains one node and all parameters to fit
##'     are in \code{ldata}, then that node will be replicated and
##'     each of the replicated nodes represent one particle in the
##'     trajectory (see \sQuote{Examples}). On the other hand if the
##'     model contains multiple nodes or the parameters to fit are
##'     contained in \code{gdata}, then the trajectory in the
##'     \code{result} argument represents one particle.
##' @template tolerance-param
##' @param ... Further arguments to be passed to \code{fn}.
##' @template verbose-param
##' @return A \code{SimInf_abc} object.
##' @references
##'
##' \Toni2009
##' @export
##' @importFrom stats cov
##' @example man/examples/abc.R
setGeneric(
    "abc",
    signature = "model",
    function(model, priors, ngen, npart, fn, tolerance, ...,
             verbose = getOption("verbose", FALSE)) {
        standardGeneric("abc")
    }
)

##' @rdname abc
##' @export
setMethod(
    "abc",
    signature(model = "SimInf_model"),
    function(model, priors, ngen, npart, fn, tolerance, ..., verbose) {
        check_integer_arg(npart)
        npart <- as.integer(npart)
        if (length(npart) != 1L || npart <= 1L)
            stop("'npart' must be an integer > 1.", call. = FALSE)

        ## Match the 'priors' to parameters in 'ldata' or 'gdata'.
        priors <- parse_priors(priors)
        pars <- match_priors(model, priors)

        object <- new("SimInf_abc", model = model, priors = priors,
                      target = pars$target, pars = pars$pars, npart = npart,
                      nprop = integer(), fn = fn, x = list(),
                      epsilon = matrix(numeric(0), ncol = 0, nrow = 0),
                      w = list(), distance = list(), ess = numeric())

        continue(object, ngen = ngen, tolerance = tolerance, ...,
                 verbose = verbose)
    }
)

##' Run more generations of ABC SMC
##'
##' @param object The \code{SimInf_abc} to continue from.
##' @param ngen The number of generations of ABC-SMC to run.
##' @template tolerance-param
##' @param ... Further arguments to be passed to
##'     \code{SimInf_abc@@fn}.
##' @template verbose-param
##' @return A \code{SimInf_abc} object.
##' @export
setGeneric(
    "continue",
    signature = "object",
    function(object, ngen = 1, tolerance, ...,
             verbose = getOption("verbose", FALSE)) {
        standardGeneric("continue")
    }
)

##' @rdname continue
##' @export
setMethod(
    "continue",
    signature(object = "SimInf_abc"),
    function(object, ngen = 1, tolerance, ...,
             verbose = getOption("verbose", FALSE)) {
        check_integer_arg(ngen)
        ngen <- as.integer(ngen)
        if (length(ngen) != 1L || ngen < 1L)
            stop("'ngen' must be an integer >= 1.", call. = FALSE)

        tolerance <- abc_tolerance(tolerance, object@tolerance, ngen)
        if (ncol(object@tolerance) == 0)
            dim(object@tolerance) <- c(nrow(tolerance), 0)

        abc_fn <- switch(object@target,
                         "gdata" = abc_gdata,
                         "ldata" = abc_ldata,
                         stop("Unknown target: ", object@target,
                              call. = FALSE))

        ## Setup a population of particles (x) and weights (w).
        x <- NULL
        if (length(object@x))
            x <- object@x[[length(object@x)]]
        w <- NULL
        if (length(object@w))
            w <- object@w[[length(object@w)]]

        ## Append new generations to object
        generations <- seq(length(object@x) + 1, length(object@x) + ngen)
        for (generation in generations) {
            tmp <- abc_fn(object@model, object@pars, object@priors,
                          object@npart, object@fn, generation,
                          epsilon, x, w, verbose, ...)

            ## Move the population of particles to the next
            ## generation.
            x <- tmp$x
            object@x[[length(object@x) + 1]] <- x
            w <- tmp$w
            object@w[[length(object@w) + 1]] <- w
            object@tolerance <- cbind(object@tolerance, tmp$tolerance)
            object@ess[length(object@ess) + 1] <- 1 / sum(w^2)
            object@nprop[length(object@nprop) + 1] <- tmp$nprop
        }

        object
    }
)

##' @rdname run
##' @include run.R
##' @export
setMethod(
    "run",
    signature(model = "SimInf_abc"),
    function(model, ...) {
        ## Sample a particle to use for the parameters from the last
        ## generation.
        generation <- length(model@x)
        particle <- sample.int(ncol(model@x[[generation]]), 1)

        ## Apply the particle to the model.
        for (i in seq_len(nrow(model@x[[generation]]))) {
            parameter <- model@pars[i]
            value <- model@x[[generation]][i, particle]
            if (identical(model@target, "gdata")) {
                model@model@gdata[parameter] <- value
            } else {
                model@model@ldata[parameter, 1] <- value
            }
        }

        ## Run the model using the particle.
        run(model@model, ...)
    }
)
