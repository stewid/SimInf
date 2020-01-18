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

##' Class \code{"SimInf_abc_smc"}
##'
##' @section Slots:
##' \describe{
##'   \item{model}{
##'     FIXME.
##'   }
##'   \item{priors}{
##'     FIXME.
##'   }
##'   \item{target}{
##'     FIXME.
##'   }
##'   \item{i}{
##'     FIXME.
##'   }
##'   \item{npart}{
##'     FIXME.
##'   }
##'   \item{fn}{
##'     FIXME.
##'   }
##'   \item{x}{
##'     FIXME.
##'   }
##'   \item{w}{
##'     FIXME.
##'   }
##' }
##' @export
setClass("SimInf_abc_smc",
         slots = c(model  = "SimInf_model",
                   priors = "data.frame",
                   target = "character",
                   i      = "integer",
                   npart  = "integer",
                   fn     = "function",
                   x      = "list",
                   w      = "list"))

##' Display the ABC posterior distribution
##'
##' @param x The \code{SimInf_abc_smc} object to plot.
##' @param y The generation to plot. The default is to display the
##'     last generation.
##' @param ... Additional arguments affecting the plot.
##' @aliases plot,SimInf_abc_smc-method
##' @importFrom graphics contour
##' @importFrom graphics lines
##' @importFrom MASS bandwidth.nrd
##' @importFrom MASS kde2d
##' @importFrom stats density
##' @export
setMethod("plot",
          signature(x = "SimInf_abc_smc"),
          function(x, y, ...) {
              if (missing(y))
                  y <- length(x@x)
              pairs(t(x@x[[y]]),
                    diag.panel = function(x, ...) {
                        usr <- par("usr")
                        on.exit(par(usr))
                        par(usr = c(usr[1:2], 0, 1.5))
                        d <- density(x, bw = "SJ-ste")
                        d$y <- d$y / max(d$y)
                        lines(d, ...)
                    },
                    lower.panel = function(x, y, ...) {
                        h <- c(bandwidth.nrd(x), bandwidth.nrd(y))
                        d <- kde2d(x, y, h = h, n = 100)
                        contour(d, add = TRUE, drawlabels = FALSE, ...)
                    }, ...)
          }
)

summary_particles <- function(x, i) {
    str <- sprintf("Particles in generation %i:", i)
    cat(sprintf("\n%s\n", str))
    cat(sprintf("%s\n", paste0(rep("-", nchar(str)), collapse = "")))
    summary_matrix(x[[i]])
}

##' Print summary of a \code{SimInf_abc_smc} object
##'
##' @param object The \code{SimInf_abc_smc} object.
##' @return \code{invisible(object)}.
##' @export
##' @importFrom methods show
setMethod("show",
          signature(object = "SimInf_abc_smc"),
          function(object) {
              cat(sprintf("Number of particles per generation: %i\n",
                          object@npart))
              cat(sprintf("Number of generations: %i\n",
                          length(object@x)))

              if (length(object@x)) {
                  summary_particles(object@x, length(object@x))
              }

              invisible(object)
          }
)

##' Detailed summary of a \code{SimInf_abc_smc} object
##'
##' @param object The \code{SimInf_abc_smc} object
##' @param ... Additional arguments affecting the summary produced.
##' @return None (invisible 'NULL').
##' @include SimInf_model.R
##' @export
setMethod("summary",
          signature(object = "SimInf_abc_smc"),
          function(object, ...) {
              cat(sprintf("Number of particles per generation: %i\n",
                          object@npart))
              cat(sprintf("Number of generations: %i\n",
                          length(object@x)))

              for (i in seq_len(length(object@x))) {
                  summary_particles(object@x, i)
              }

              invisible(NULL)
          }
)

##' Check model before running ABC-SMC
##'
##' Raise an error if the model argument is not ok.
##' @param model the model to check.
##' @return invisible(NULL)
##' @noRd
check_model_for_abc_smc <- function(model) {
    check_model_argument(model)

    if (!identical(Nn(model), 1L))
        stop("The 'model' must contain one node.", call. = FALSE)
    if (length(model@events@event) > 0)
        stop("The 'model' cannot contain any events.", call. = FALSE)

    invisible(NULL)
}

##' Generate replicates of first node in model.
##'
##' Replicate the node specific matrices 'u0', 'v0' and 'ldata' in the
##' first node.
##' @param model the model to replicate.
##' @param n the number of replicates.
##' @return A modified model object
##' @noRd
replicate_first_node <- function(model, n) {
    if (dim(model@u0)[1] > 0)
        model@u0 <- model@u0[, rep(1, n), drop = FALSE]
    if (dim(model@v0)[1] > 0)
        model@v0 <- model@v0[, rep(1, n), drop = FALSE]
    if (dim(model@ldata)[1] > 0)
        model@ldata <- model@ldata[, rep(1, n), drop = FALSE]
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

##' @importFrom utils setTxtProgressBar
##' @importFrom utils txtProgressBar
##' @noRd
abc_smc_gdata <- function(model, i, priors, npart, fn,
                          generation, x, w, verbose, ...) {
    if (isTRUE(verbose)) {
        cat("\nGeneration", generation, "...\n")
        pb <- txtProgressBar(min = 0, max = npart, style = 3)
        t0 <- proc.time()
    }

    xx <- NULL
    ancestor <- NULL
    tot_proposals <- 0
    sigma <- proposal_covariance(x)

    while (n_particles(xx) < npart) {
        proposals <- .Call(SimInf_abc_smc_proposals,
                           priors$parameter, priors$distribution,
                           priors$p1, priors$p2, 1L, x, w, sigma)
        for (j in seq_len(nrow(proposals))) {
            model@gdata[i[j]] <- proposals[j, 1]
        }

        result <- fn(run(model), generation, ...)
        stopifnot(is.logical(result), length(result) == 1L)
        tot_proposals <- tot_proposals + 1L
        if (isTRUE(result)) {
            ## Collect accepted particle
            xx <- cbind(xx, as.matrix(model@gdata)[i, 1, drop = FALSE])
            ancestor <- c(ancestor, attr(proposals, "ancestor")[1])
        }

        ## Report progress.
        if (isTRUE(verbose)) {
            setTxtProgressBar(pb, n_particles(xx))
        }
    }

    ## Calculate weights.
    ww <- .Call(SimInf_abc_smc_weights, priors$distribution,
                priors$p1, priors$p2, x[, ancestor], xx, w, sigma)

    ## Report progress.
    if (isTRUE(verbose)) {
        t1 <- proc.time()
        cat(sprintf("\n\n  accrate = %.2e, ESS = %.2e time = %.2f secs\n\n",
                    npart / tot_proposals, 1 / sum(ww^2), (t1 - t0)[3]))
        summary_matrix(xx)
    }

    list(x = xx, w = ww)
}

##' @importFrom utils setTxtProgressBar
##' @importFrom utils txtProgressBar
##' @noRd
abc_smc_ldata <- function(model, i, priors, npart, fn,
                          generation, x, w, verbose, ...) {
    ## Let each node represents one particle. Replicate the first node
    ## to run many particles simultanously. Start with 10 x 'npart'
    ## and then increase the number adaptively based on the acceptance
    ## rate.
    n <- as.integer(10 * npart)
    model <- replicate_first_node(model, n)

    if (isTRUE(verbose)) {
        cat("\nGeneration", generation, "...\n")
        pb <- txtProgressBar(min = 0, max = npart, style = 3)
        t0 <- proc.time()
    }

    xx <- NULL
    ancestor <- NULL
    tot_proposals <- 0
    sigma <- proposal_covariance(x)

    while (n_particles(xx) < npart) {
        if (all(n < 1e5, tot_proposals > 2 * n)) {
            ## Increase the number of particles that is simulated in
            ## each trajectory.
            n <- min(1e5L, n * 2L)
            model <- replicate_first_node(model, n)
        }

        proposals <- .Call(SimInf_abc_smc_proposals,
                           priors$parameter, priors$distribution,
                           priors$p1, priors$p2, n, x, w, sigma)
        for (j in seq_len(nrow(proposals))) {
            model@ldata[i[j], ] <- proposals[j, ]
        }

        result <- fn(run(model), generation, ...)
        stopifnot(is.logical(result), length(result) == n)

        ## Collect accepted particles making sure not to collect more
        ## than 'npart'.
        j <- cumsum(result) + n_particles(xx)
        j <- which(j == npart)
        result <- which(result)
        if (length(j)) {
            j <- min(j)
            result <- result[result <= j]
            tot_proposals <- tot_proposals + j
            xx <- cbind(xx, model@ldata[i, result, drop = FALSE])
            ancestor <- c(ancestor, attr(proposals, "ancestor")[result])
        } else {
            tot_proposals <- tot_proposals + n
            if (length(result)) {
                xx <- cbind(xx, model@ldata[i, result, drop = FALSE])
                ancestor <- c(ancestor, attr(proposals, "ancestor")[result])
            }
        }

        ## Report progress.
        if (isTRUE(verbose)) {
            setTxtProgressBar(pb, n_particles(xx))
        }
    }

    ## Calculate weights.
    ww <- .Call(SimInf_abc_smc_weights, priors$distribution,
                priors$p1, priors$p2, x[, ancestor], xx, w, sigma)

    ## Report progress.
    if (isTRUE(verbose)) {
        t1 <- proc.time()
        cat(sprintf("\n\n  accrate = %.2e, ESS = %.2e time = %.2f secs\n\n",
                    npart / tot_proposals, 1 / sum(ww^2), (t1 - t0)[3]))
        summary_matrix(xx)
    }

    list(x = xx, w = ww)
}

##' Run ABC SMC
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
##'     the simulated trajectory and determine for each node
##'     (particle) if it should be accepted (\code{TRUE}) or rejected
##'     (\code{FALSE}). The first argument in \code{fn} is the
##'     simulated model containing one trajectory. The second argument
##'     to \code{fn} is an integer with the \code{generation}
##'     particles. The function should return a logical vector with
##'     one value for each node in the simulated model.
##' @param ... Further arguments to be passed to \code{fn}.
##' @param verbose = TRUE
##' @return A \code{SimInf_abc_smc} object.
##' @export
##' @importFrom stats cov
##' @examples
##' ## Fit a model to data from an influenza in a boarding school in
##' ## England (Anonymous. 1978. Influenza in a boarding school.
##' ## British Medical Journal 1:578.)
##' obs <- data.frame(time = 1:15,
##'                   R1 = c(0, 1, 6, 26, 73, 222, 293, 258,
##'                          236, 191, 124, 69, 26, 11, 4))
##'
##' ## The distance function to accept or reject a proposal. Each node
##' ## in the simulated trajectory (contained in the 'result' object)
##' ## represents one proposal. The 'generation' argument is the current
##' ## generation of proposals.
##' acceptFun <- function(result, generation, tol, ptol, ...) {
##'     ## Determine the tolerance for this generation.
##'     tol <- tol * (ptol)^(generation - 1)
##'
##'     ## Extract the time-series for R1 for each node as a
##'     ## data.frame.
##'     sim <- trajectory(result, "R1")
##'
##'     ## Split the 'sim' data.frame by node and calculate the sum
##'     ## of the squared distance at each time-point for each node.
##'     dist <- tapply(sim$R1, sim$node, function(sim_R1) {
##'         sum((obs$R1 - sim_R1)^2)
##'     })
##'
##'     ## Return TRUE or FALSE for each node depending on if the
##'     ## distance is less than the tolerance.
##'     dist < tol
##' }
##'
##' ## Transitions
##' transitions <- c("S -> beta*S*I/(S+I+R1+R2) -> I",
##'                  "I -> gamma1*I -> R1",
##'                  "R1 -> gamma2*R1 -> R2")
##'
##' ## Specify the compartments.
##' compartments <- c("S", "I", "R1", "R2")
##'
##' ## Create the model.
##' model <- mparse(transitions = transitions,
##'                 compartments = compartments,
##'                 ldata = data.frame(beta = 0.1, gamma1 = 0.1, gamma2 = 0.1),
##'                 u0 = data.frame(S = 762, I = 1, R1 = 0, R2 = 0),
##'                 tspan = 1:15)
##'
##' ## Fit the model parameters using ABC-SMC. The priors for the paramters
##' ## are specified in the second argument using a formula notation. Here
##' ## we use a uniform distribtion for each parameter with lower bound = 0
##' ## and upper bound = 5.
##' fit <- abc_smc(model = model,
##'                priors = c(beta~U(0, 5), gamma1~U(0, 5), gamma2~U(0, 5)),
##'                ngen = 2,
##'                npart = 50,
##'                fn = acceptFun,
##'                tol = 100000,
##'                ptol = 0.5)
##'
##' plot(fit)
##'
##' fit <- continue(fit, tol = 100000, ptol = 0.5)
##'
##' plot(fit)
abc_smc <- function(model, priors, ngen, npart, fn, ..., verbose = TRUE) {
    check_model_for_abc_smc(model)
    check_integer_arg(ngen, npart)
    ngen <- as.integer(ngen)
    npart <- as.integer(npart)

    ## Match the 'priors' to parameters in 'ldata' or 'gdata'.
    priors <- parse_priors(priors)
    i <- match(priors$parameter, rownames(model@ldata))
    if (any(is.na(i))) {
        i <- match(priors$parameter, names(model@gdata))
        if (any(is.na(i))) {
            stop("All parameters in 'priors' must be either ",
                 "in 'gdata' or 'ldata'", call. = FALSE)
        }
        target <- "gdata"
    } else {
        target <- "ldata"
    }

    object <- new("SimInf_abc_smc", model = model,
                  priors = priors, target = target,
                  i = i, npart = npart, fn = fn,
                  x = list(), w = list())

    continue(object, ngen = ngen, verbose = verbose, ...)
}

##' Run more generations of ABC SMC
##'
##' @param object The \code{SimInf_abc_smc} to continue from.
##' @param ngen The number of generations of ABC-SMC to run.
##' @param ... Further arguments to be passed to
##'     \code{SimInf_abc_smc@@fn}.
##' @param verbose = TRUE
##' @return A \code{SimInf_abc_smc} object.
##' @export
continue <- function(object, ngen = 1, ..., verbose = TRUE) {
    stopifnot(inherits(object, "SimInf_abc_smc"))
    check_integer_arg(ngen)

    abc_smc_fn <- switch(object@target,
                         "gdata" = abc_smc_gdata,
                         "ldata" = abc_smc_ldata,
                         stop("Unknown target: ", object@target,
                              call. = FALSE))

    ## Setup a population of particles (x) and weights (w)
    x <- NULL
    if (length(object@x))
        x <- object@x[[length(object@x)]]
    w <- NULL
    if (length(object@w))
        w <- object@w[[length(object@w)]]

    for (generation in seq(length(object@x) + 1, length(object@x) + ngen)) {
        tmp <- abc_smc_fn(object@model, object@i, object@priors, object@npart,
                          object@fn, generation, x, w, verbose, ...)

        ## Move the population of particles to the next generation.
        x <- tmp$x
        w <- tmp$w
        object@x[[length(object@x) + 1]] <- x
        object@w[[length(object@w) + 1]] <- w
    }

    object
}
