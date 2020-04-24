## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2020 Stefan Widgren
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
##' @slot priors A \code{data.frame} containing the four columns
##'     \code{parameter}, \code{distribution}, \code{p1} and
##'     \code{p2}. The column \code{parameter} gives the name of the
##'     parameter referred to in the model. The column
##'     \code{distribution} contains a letter indicating the prior
##'     distribution. Valid letters are 'G' (gamma), 'N' (normal) or
##'     'U' (uniform). The column \code{p1} is a numeric vector with
##'     the first hyperparameter for each prior: 'G') shape, 'N')
##'     mean, and 'U') lower bound. The column \code{p2} is a numeric
##'     vector with the second hyperparameter for each prior: 'G')
##'     rate, 'N') standard deviation, and 'U') upper bound.
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
##' @slot epsilon A numeric matrix (number of summary statistics X
##'     number of generations) where each column contains the
##'     tolerances for a generation and each row contains a sequence
##'     of gradually decreasing tolerances.
##' @slot x A list where each item is a \code{matrix} with the
##'     accepted particles in each generation. Each column is one
##'     particle.
##' @slot w A list where each item is a vector with the weights for
##'     the particles \code{x} in the corresponding generation.
##' @slot ess A numeric vector with the effective sample size (ESS) in
##'     each generation. Effective sample size is computed as
##'     \deqn{\left(\sum_{i=1}^N\!(w_{g}^{(i)})^2\right)^{-1},}{1 /
##'     (sum(w_ig^2)),} where \eqn{w_{g}^{(i)}}{w_ig} is the
##'     normalized weight of particle \eqn{i} in generation \eqn{g}.
##' @seealso \code{\link{abc}} and \code{\link{continue}}.
##' @export
setClass("SimInf_abc",
         slots = c(model   = "SimInf_model",
                   priors  = "data.frame",
                   target  = "character",
                   pars    = "integer",
                   npart   = "integer",
                   nprop   = "integer",
                   fn      = "function",
                   epsilon = "matrix",
                   x       = "list",
                   w       = "list",
                   ess     = "numeric"))

setAs(from = "SimInf_abc",
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
setMethod("show",
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
setMethod("summary",
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
##' @param n_events the number of of events for the first node in the
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

abc_progress <- function(t0, t1, x, w, npart, nprop) {
    t1 <- proc.time()
    cat(sprintf("\n\n  accrate = %.2e, ESS = %.2e time = %.2f secs\n\n",
                npart / nprop, 1 / sum(w^2), (t1 - t0)[3]))
    summary_matrix(x)
}

##' Check that the returned result from the abc distance function is
##' valid.
##' @noRd
check_abc_accept <- function(result, n, old_epsilon, epsilon) {
    if (!is.list(result)) {
        stop("The result from the ABC distance function must be a 'list'.",
             call. = FALSE)
    }

    if (!is.logical(result$accept)) {
        stop("The accepted/rejected vector must be of type 'logical'.",
             call. = FALSE)
    }

    if (!identical(length(result$accept), n)) {
        stop("Invalid length of the  accepted/rejected vector.",
             call. = FALSE)
    }

    if (!all(is.vector(result$epsilon, "numeric"),
             length(result$epsilon) > 0,
             all(result$epsilon >= 0))) {
        stop("'epsilon' must be a numeric vector with non-negative values.",
             call. = FALSE)
    }

    if (is.null(epsilon)) {
        if (!is.null(old_epsilon)) {
            if (!all(result$epsilon < old_epsilon)) {
                stop("'epsilon' must decrease in each generation.",
                     call. = FALSE)
            }
        }

        return(TRUE)
    }

    if (!identical(epsilon, result$epsilon)) {
        stop("'epsilon' must be fixed within each generation.",
             call. = FALSE)
    }

    FALSE
}

##' Return result for ABC whether the particles match the data
##'
##' @param accept A logical vector with one value for each particle in
##'     the simulated model.
##' @param epsilon A numeric vector with the tolerance for the
##'     generation.
##' @return \code{list(accept = accept, epsilon = epsilon)}
##' @export
abc_accept <- function(accept, epsilon) {
    list(accept = accept, epsilon = epsilon)
}

##' @importFrom utils setTxtProgressBar
##' @importFrom utils txtProgressBar
##' @noRd
abc_gdata <- function(model, pars, priors, npart, fn, generation,
                      old_epsilon, x, w, verbose, ...) {
    if (isTRUE(verbose)) {
        cat("\nGeneration", generation, "...\n")
        pb <- txtProgressBar(min = 0, max = npart, style = 3)
        t0 <- proc.time()
    }

    xx <- NULL
    ancestor <- NULL
    epsilon <- NULL
    nprop <- 0L
    sigma <- proposal_covariance(x)

    while (n_particles(xx) < npart) {
        proposals <- .Call(SimInf_abc_proposals, priors$parameter,
                           priors$distribution, priors$p1, priors$p2,
                           1L, x, w, sigma)
        for (i in seq_len(nrow(proposals))) {
            model@gdata[pars[i]] <- proposals[i, 1]
        }

        result <- fn(run(model), generation, ...)
        if (check_abc_accept(result, 1L, old_epsilon, epsilon))
            epsilon <- result$epsilon
        nprop <- nprop + 1L
        if (isTRUE(result$accept)) {
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

    list(x = xx, w = ww, nprop = nprop, epsilon = epsilon)
}

##' @importFrom utils setTxtProgressBar
##' @importFrom utils txtProgressBar
##' @noRd
abc_ldata <- function(model, pars, priors, npart, fn, generation,
                      old_epsilon, x, w, verbose, ...) {
    ## Let each node represents one particle. Replicate the first node
    ## to run many particles simultanously. Start with 10 x 'npart'
    ## and then increase the number adaptively based on the acceptance
    ## rate.
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
    epsilon <- NULL
    nprop <- 0L
    sigma <- proposal_covariance(x)

    while (n_particles(xx) < npart) {
        if (all(n < 1e5L, nprop > 2L * n)) {
            ## Increase the number of particles that is simulated in
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

        result <- fn(run(model), generation, ...)
        if (check_abc_accept(result, n, old_epsilon, epsilon))
            epsilon <- result$epsilon

        ## Collect accepted particles making sure not to collect more
        ## than 'npart'.
        i <- cumsum(result$accept) + n_particles(xx)
        i <- which(i == npart)
        j <- which(result$accept)
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

    list(x = xx, w = ww, nprop = nprop, epsilon = epsilon)
}

##' Approximate Bayesian computation
##'
##' @param model The model to generate data from.
##' @param priors The priors for the parameters to fit. Each prior is
##'     specified with a formula notation, for example, \code{beta ~
##'     U(0, 1)} to specify that beta is uniformly distributed between
##'     0 and 1. Use \code{c()} to provide more than one prior, for
##'     example, \code{c(beta ~ U(0, 1), gamma ~ N(10, 1)}. Gamma
##'     \code{G}, normal \code{N} and uniform \code{U} distributions
##'     are supported.
##' @param ngen The number of generations of ABC-SMC to run.
##' @param npart An integer specifying the number of particles.
##' @param fn A function for calculating the summary statistics for
##'     the simulated trajectory and determine for each particle if it
##'     should be accepted (\code{TRUE}) or rejected (\code{FALSE}).
##'     The first argument in \code{fn} is the simulated model
##'     containing one trajectory.  The second argument to \code{fn}
##'     is an integer with the \code{generation} of the particles.
##'     The function should return a logical vector with one value for
##'     each particle in the simulated model.
##' @param ... Further arguments to be passed to \code{fn}.
##' @template verbose-param
##' @return A \code{SimInf_abc} object.
##' @export
##' @importFrom stats cov
##' @example man/examples/abc.R
abc <- function(model, priors, ngen, npart, fn, ...,
                verbose = getOption("verbose", FALSE)) {
    check_model_argument(model)
    check_integer_arg(npart)
    npart <- as.integer(npart)
    if (length(npart) != 1L || npart <= 1L)
        stop("'npart' must be an integer > 1.", call. = FALSE)

    ## Match the 'priors' to parameters in 'ldata' or 'gdata'.
    priors <- parse_priors(priors)
    pars <- match(priors$parameter, rownames(model@ldata))
    if (any(is.na(pars))) {
        pars <- match(priors$parameter, names(model@gdata))
        if (any(is.na(pars))) {
            stop("All parameters in 'priors' must be either ",
                 "in 'gdata' or 'ldata'.", call. = FALSE)
        }
        target <- "gdata"
    } else {
        if (!identical(Nn(model), 1L))
            stop("The 'model' must contain one node.", call. = FALSE)
        target <- "ldata"
    }

    object <- new("SimInf_abc", model = model, priors = priors,
                  target = target, pars = pars, npart = npart,
                  nprop = integer(), fn = fn, x = list(),
                  epsilon = matrix(numeric(0), ncol = 0, nrow = 0),
                  w = list(), ess = numeric())

    continue(object, ngen = ngen, verbose = verbose, ...)
}

##' Run more generations of ABC SMC
##'
##' @param object The \code{SimInf_abc} to continue from.
##' @param ngen The number of generations of ABC-SMC to run.
##' @param ... Further arguments to be passed to
##'     \code{SimInf_abc@@fn}.
##' @template verbose-param
##' @return A \code{SimInf_abc} object.
##' @export
continue <- function(object, ngen = 1, ...,
                     verbose = getOption("verbose", FALSE)) {
    stopifnot(inherits(object, "SimInf_abc"))
    check_integer_arg(ngen)
    ngen <- as.integer(ngen)
    if (length(ngen) != 1L || ngen < 1L)
        stop("'ngen' must be an integer >= 1.", call. = FALSE)

    abc_fn <- switch(object@target,
                     "gdata" = abc_gdata,
                     "ldata" = abc_ldata,
                     stop("Unknown target: ", object@target,
                          call. = FALSE))

    ## Setup a population of particles (x), weights (w) and epsilon.
    x <- NULL
    if (length(object@x))
        x <- object@x[[length(object@x)]]
    w <- NULL
    if (length(object@w))
        w <- object@w[[length(object@w)]]
    epsilon <- NULL
    if (ncol(object@epsilon))
        epsilon <- object@epsilon[, ncol(object@epsilon)]

    ## Append new generations to object
    generations <- seq(length(object@x) + 1, length(object@x) + ngen)
    for (generation in generations) {
        tmp <- abc_fn(object@model, object@pars, object@priors,
                      object@npart, object@fn, generation,
                      epsilon, x, w, verbose, ...)

        ## Move the population of particles to the next generation.
        x <- tmp$x
        object@x[[length(object@x) + 1]] <- x
        w <- tmp$w
        object@w[[length(object@w) + 1]] <- w
        epsilon <- tmp$epsilon
        if (ncol(object@epsilon) == 0)
            dim(object@epsilon) <- c(length(tmp$epsilon), 0)
        object@epsilon <- cbind(object@epsilon, epsilon)
        object@ess[length(object@ess) + 1] <- 1 / sum(w^2)
        object@nprop[length(object@nprop) + 1] <- tmp$nprop
    }

    object
}
