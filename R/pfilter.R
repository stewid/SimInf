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

##' Class \code{"SimInf_pfilter"}
##'
##' @slot model A \code{SimInf_model} object with one filtered
##'     trajectory attached.
##' @slot npart An integer with the number of particles that was used
##'     at each timestep.
##' @slot loglik The estimated log likelihood.
##' @slot ess A numeric vector with the effective sample size (ESS).
##'     The effective sample size is computed as
##'     \deqn{\left(\sum_{i=1}^N\!(w_{t}^{i})^2\right)^{-1},}{1 /
##'     (sum(w_it^2)),} where \eqn{w_{t}^{i}}{w_it} is the normalized
##'     weight of particle \eqn{i} at time \eqn{t}.
##' @export
setClass(
    "SimInf_pfilter",
    slots = c(model  = "SimInf_model",
              npart  = "integer",
              loglik = "numeric",
              ess    = "numeric")
)

##' Brief summary of a \code{SimInf_pfilter} object
##'
##' @param object The \code{SimInf_pfilter} object.
##' @return \code{invisible(object)}.
##' @export
setMethod(
    "show",
    signature(object = "SimInf_pfilter"),
    function(object) {
        cat(sprintf("Number of particles: %i\n", object@npart))
        cat(sprintf("Log-likelihood: %f\n", object@loglik))

        invisible(object)
    }
)

##' Detailed summary of a \code{SimInf_pfilter} object
##'
##' @param object The \code{SimInf_pfilter} object.
##' @param ... Unused additional arguments.
##' @return \code{invisible(NULL)}.
##' @export
setMethod(
    "summary",
    signature(object = "SimInf_pfilter"),
    function(object, ...) {
        cat("Particle filter\n")
        cat("---------------\n")
        cat(sprintf("Number of particles: %i\n", object@npart))
        cat(sprintf("Log-likelihood: %f\n", object@loglik))
        cat(sprintf("Model: %s\n", as.character(class(object@model))))
        cat(sprintf("Number of nodes: %i\n", n_nodes(object@model)))
        show(object@model@events)

        summary_transitions(object@model)

        if (length(object@model@gdata))
            summary_gdata(object@model)
        if (ncol(object@model@ldata))
            summary_data_matrix(object@model@ldata, "Local data")
        summary_output_matrix(object@model, "Continuous state variables",
                              rownames(object@model@v0))
        summary_output_matrix(object@model, "Compartments",
                              rownames(object@model@S))

        invisible(NULL)
    }
)

##' Split data into intervals.
##'
##' @return A list of data.frames where each list item is the data for
##'     one interval in tspan.
##' @noRd
pfilter_data <- function(model, data) {
    if (!is.data.frame(data))
        stop("'data' must be a data.frame.", call. = FALSE)

    if (!("time" %in% names(data)))
        stop("Missing 'time' column in data.", call. = FALSE)

    if (!is.null(names(model@tspan))) {
        data$time <- julian(x = as.Date(data$time),
                            origin = as.Date(names(model@tspan)[1]))
        data$time <- as.numeric(data$time + model@tspan[1])
    }

    check_integer_arg(data$time)
    data <- data[order(data$time), seq_len(ncol(data)), drop = FALSE]

    if (any(length(data$time) < 1,
            anyNA(data$time),
            any(diff(unique(data$time)) <= 0))) {
        stop("'time' column in data must be an increasing vector.",
             call. = FALSE)
    }

    if (data$time[1] < model@tspan[1])
        stop("data$time[1] must be >= tspan[1].", call. = TRUE)

    lapply(split(data, data$time), function(x) {
        rownames(x) <- NULL
        x
    })
}

##' Split tspan into intervals.
##'
##' @return A two-column matrix where each row specifies the tspan to
##'     use when running the model from time[i] to time[i+1]. The
##'     first column is \code{NA} if the interval is one time-unit.
##' @noRd
pfilter_tspan <- function(model, data) {
    time <- sapply(data, function(x) {
        x$time[1]
    })

    do.call("rbind", lapply(seq_len(length(time)), function(i) {
        if (i == 1) {
            if (model@tspan[1] < time[1])
                return(as.numeric(c(model@tspan[1], time[1])))
            return(as.numeric(c(NA_real_, time[i])))
        }

        if (diff(as.integer(time[c(i - 1L, i)])) > 1)
            return(as.numeric(time[c(i - 1L, i)]))
        as.numeric(c(NA_real_, time[i]))
    }))
}

##' Split scheduled events into the intervals in tspan.
##'
##' @param events The scheduled events to split.
##' @param time_end Endpoint time in each tpsan interval.
##' @return A list with the scheduled events to use in each interval.
##' @noRd
pfilter_events <- function(events, time_end) {
    if (length(events@event) == 0)
        return(NULL)

    m <- .Call(SimInf_split_events, events@time, as.integer(time_end))
    lapply(seq_len(nrow(m)), function(i) {
        j <- seq(from = m[i, 1], by = 1, length.out = m[i, 2])

        methods::new("SimInf_events",
                     E          = events@E,
                     N          = events@N,
                     event      = events@event[j],
                     time       = events@time[j],
                     node       = events@node[j],
                     dest       = events@dest[j],
                     n          = events@n[j],
                     proportion = events@proportion[j],
                     select     = events@select[j],
                     shift      = events@shift[j])
    })
}

pfilter_obs_process <- function(model, obs_process, data, npart) {
    if (is.function(obs_process))
        return(match.fun(obs_process))

    if (n_nodes(model) > 1) {
        stop("The observation process must be a function ",
             "for a model with multiple nodes.",
             call. = FALSE)
    }

    if (any(as.integer(sapply(data, nrow)) > 1L)) {
        stop("The observation process must be a function ",
             "when data contains multiple rows for a ",
             "time-point.", call. = FALSE)
    }

    if (!methods::is(obs_process, "formula")) {
        stop("'obs_process' must be either a formula or a function.",
             call. = FALSE)
    }

    obs_process <- parse_distribution(obs_process)

    ## Match the parameter on the lhs of the observation process to a
    ## column in the data data.frame.
    data_columns <- colnames(data[[1]])
    par_i <- match(obs_process$parameter, data_columns)
    par <- data_columns[par_i]
    if (!isTRUE(par %in% data_columns)) {
        stop("Unable to match the parameter on the lhs to a column in 'data'.",
             call. = FALSE)
    }

    ## Match the symbols on the rhs of the observation process to
    ## compartments in U or V.
    symbols <- unlist(obs_process$symbols)
    u <- lapply(match(symbols, rownames(model@S), nomatch = 0), function(i) {
        list(slot = "U",
             name = rownames(model@S)[i],
             i = i)
    })

    symbols <- setdiff(symbols, sapply(u, function(x) x$name))
    v <- lapply(match(symbols, rownames(model@v0), nomatch = 0), function(i) {
        list(slot = "V",
             name = rownames(model@v0)[i],
             i = i)
    })

    symbols <- setdiff(symbols, sapply(v, function(x) x$name))
    if (length(symbols) > 0) {
        stop("Non-existing compartment(s) in model: ",
             paste0("'", symbols, "'", collapse = ", "),
             ".", call. = FALSE)
    }

    expr <- switch(obs_process$distribution,
                   binomial = {
                       paste0("stats::dbinom(x = ",
                              obs_process$parameter,
                              ", size = ",
                              obs_process$p1,
                              ", prob = ",
                              obs_process$p2,
                              ", log = TRUE)")
                   },
                   poisson = {
                       paste0("stats::dpois(x = ",
                              obs_process$parameter,
                              ", lambda = ",
                              obs_process$p1,
                              ", log = TRUE)")
                   },
                   uniform = {
                       paste0("stats::dunif(x = ",
                              obs_process$parameter,
                              ", min = ",
                              obs_process$p1,
                              ", max = ",
                              obs_process$p2,
                              ", log = TRUE)")
                   },
                   stop("Unknown distribution: '",
                        obs_process$distribution, "'.",
                        call. = FALSE)
                   )

    list(slots = c(u, v), expr = expr, par = par, par_i = par_i)
}

##' Run a particle filter on a model.
##' @noRd
pfilter_internal <- function(model,
                             events,
                             obs_process,
                             data,
                             npart,
                             tspan,
                             init_model) {
    Nd <- nrow(model@v0)
    Ntspan <- nrow(tspan)
    ess <- numeric(Ntspan)
    loglik <- 0
    w <- numeric(npart)

    ## Create a matrix to keep track of the states in the compartments
    ## (including the initial state) in every particle.
    U <- matrix(data = NA_integer_,
                nrow = npart * n_nodes(model) * n_compartments(model),
                ncol = Ntspan + 1L)

    ## Create a matrix to keep track of the continuous states
    ## (including the initial state) in every particle.
    V <- matrix(data = NA_real_,
                nrow = npart * n_nodes(model) * Nd,
                ncol = Ntspan + 1L)

    ## Initialise the model.
    if (is.function(init_model))
        model <- init_model(model)
    U[, 1L] <- rep(as.integer(model@u0), npart)
    V[, 1L] <- rep(as.numeric(model@v0), npart)

    a <- matrix(data = NA_integer_,
                nrow = npart,
                ncol = Ntspan + 1L)
    a[, 1L] <- seq_len(npart)

    ## Loop over time series.
    m <- model
    m@replicates <- npart
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
            rep((a[, i] - 1L) * n_nodes(model) * n_compartments(model),
                each = n_nodes(model) * n_compartments(model))
        m@u0 <- matrix(
            data = U[u_i, i],
            nrow = n_compartments(model),
            dimnames = dimnames(m@u0))

        v_i <- seq_len(n_nodes(model) * Nd) +
            rep((a[, i] - 1L) * n_nodes(model) * Nd,
                each = n_nodes(model) * Nd)
        m@v0 <- matrix(
            data = V[v_i, i],
            nrow = nrow(m@v0),
            dimnames = dimnames(m@v0))

        ## Propagate the model.
        x <- run(m)
        if (length(x@tspan) > 1L) {
            ## Extract the final state in each replicate. The first
            ## replicate is in column 1 and 2, the second replicate in
            ## column 3 and 4, et cetera.
            x@tspan <- x@tspan[2L]
            j <- seq.int(from = 2L, by = 2L, length.out = npart)
            x@U <- x@U[, j, drop = FALSE]
            x@V <- x@V[, j, drop = FALSE]
        }

        ## Save states.
        U[, i + 1L] <- x@U
        V[, i + 1L] <- x@V

        ## Set the weight for the particle.
        if (is.function(obs_process)) {
            w <- obs_process(x, data[[i]])
        } else {
            e <- new.env(parent = baseenv())

            assign(x = obs_process$par,
                   value = data[[i]][, obs_process$par_i],
                   pos = e)

            for (j in seq_len(length(obs_process$slots))) {
                s <- obs_process$slots[[j]]$slot
                v <- slot(x, s)[obs_process$slots[[j]]$i, ]
                n <- obs_process$slots[[j]]$name
                assign(x = n, value = v, pos = e)
            }

            expr <- obs_process$expr
            e$expr <- expr
            w <- evalq(eval(parse(text = expr)), envir = e)
        }

        if (!all(all(is.finite(w)),
                 identical(length(w), npart),
                 is.vector(w, "numeric"))) {
            stop("Invalid observation process.", call. = FALSE)
        }

        max_w <- max(w)
        w <- exp(w - max_w)
        sum_w <- sum(w)
        loglik <- loglik + max_w + log(sum_w) - log(npart)
        w <- w / sum_w
        ess[i] <- 1 / sum(w^2)

        ## Resampling
        a[, i + 1L] <- .Call(SimInf_systematic_resampling, w)
    }

    ## Sample a trajectory.
    model@U <- matrix(data = NA_integer_,
                      nrow = n_nodes(model) * n_compartments(model),
                      ncol = Ntspan)
    model@V <- matrix(data = NA_real_,
                      nrow = n_nodes(model) * Nd,
                      ncol = Ntspan)

    i <- sample.int(npart, 1)
    for (j in rev(seq_len(Ntspan))) {
        u_i <- seq.int(
            from = (i - 1L) * n_nodes(model) * n_compartments(model) + 1L,
            length.out = n_nodes(model) * n_compartments(model))
        model@U[, j] <- U[u_i, j + 1L, drop = FALSE]

        v_i <- seq.int(
            from = (i - 1L) * n_nodes(model) * Nd + 1L,
            length.out = n_nodes(model) * Nd)
        model@V[, j] <- V[v_i, j + 1L, drop = FALSE]

        i <- a[i, j]
    }

    u_i <- seq.int(
        from = (i - 1L) * n_nodes(model) * n_compartments(model) + 1L,
        length.out = n_nodes(model) * n_compartments(model))
    model@u0 <- matrix(
        data = U[u_i, 1L],
        nrow = n_compartments(model),
        dimnames = dimnames(model@u0))

    v_i <- seq.int(
        from = (i - 1L) * n_nodes(model) * Nd + 1L,
        length.out = n_nodes(model) * Nd)
    m@v0 <- matrix(
        data = V[v_i, 1L],
        nrow = nrow(model@v0),
        dimnames = dimnames(model@v0))

    methods::new("SimInf_pfilter",
                 model = model,
                 npart = npart,
                 loglik = loglik,
                 ess = ess)
}

##' Bootstrap particle filter
##'
##' The bootstrap filtering algorithm. Systematic resampling is
##' performed at each observation.
##'
##' @param model The \code{SimInf_model} object to simulate data from.
##' @template obs_process-param
##' @template data-param
##' @template npart-param
##' @template init_model-param
##' @return A \code{SimInf_pfilter} object.
##' @references
##'
##' \Gordon1993
##' @export
##' @example man/examples/pfilter.R
setGeneric(
    "pfilter",
    signature = "model",
    function(model,
             obs_process,
             data,
             npart,
             init_model = NULL) {
        standardGeneric("pfilter")
    }
)

##' @rdname pfilter
##' @export
setMethod(
    "pfilter",
    signature(model = "SimInf_model"),
    function(model,
             obs_process,
             data,
             npart,
             init_model) {
        if (isFALSE(identical(dim(model@U_sparse), c(0L, 0L))) ||
            isFALSE(identical(dim(model@V_sparse), c(0L, 0L)))) {
            stop("'pfilter' cannot run a model with a sparse result matrix.",
                 call. = FALSE)
        }

        npart <- check_n_particles(npart)
        data <- pfilter_data(model, data)
        tspan <- pfilter_tspan(model, data)
        model@tspan <- tspan[, 2]
        events <- pfilter_events(model@events, tspan[, 2] - 1)
        obs_process <- pfilter_obs_process(model, obs_process, data, npart)
        if (!is.null(init_model))
            init_model <- match.fun(init_model)

        pfilter_internal(model,
                         events,
                         obs_process,
                         data,
                         npart,
                         tspan,
                         init_model)
   }
)

##' Diagnostic plot of a particle filter object
##'
##' @param x The \code{SimInf_pfilter} object to plot.
##' @param y If y is \code{NULL} or missing (default), the filtered
##'     trajectory (top) and the effective sample size (bottom) are
##'     displayed. If \code{y} is a character vector or a formula, the
##'     plot function for a \code{SimInf_model} object is called with
##'     the filtered trajectory, see
##'     \code{\link{plot,SimInf_model-method}} for more details about
##'     the specification a plot.
##' @param ... Other graphical parameters that are passed on to the
##'     plot function.
##' @aliases plot,SimInf_pfilter-method
##' @export
setMethod(
    "plot",
    signature(x = "SimInf_pfilter"),
    function(x, y, ...) {
        if (missing(y)) {
            savepar <- par(mfrow = c(2, 1))
            on.exit(par(savepar), add = TRUE)

            ## Settings for the x-axis
            if (is.null(names(x@model@tspan))) {
                xx <- x@model@tspan
                xlab <- "Time"
            } else {
                xx <- as.Date(names(x@model@tspan))
                xlab <- "Date"
            }

            ## Plot the sampled trajectory.
            plot(x@model, ...)

            ## Plot the effective sample size.
            plot(xx, x@ess, xlab = xlab, ylab = "ESS",
                 ylim = c(0, x@npart), frame.plot = FALSE, type = "l")
        } else {
            ## Plot the sampled trajectory.
            plot(x@model, y, ...)
        }

        invisible(NULL)
    }
)

##' Log likelihood
##'
##' Extract the estimated log likelihood from a \code{SimInf_pfilter}
##' object.
##'
##' @param object The \code{SimInf_pfilter} object.
##' @return the estimated log likelihood.
##' @export
setMethod(
    "logLik",
    signature(object = "SimInf_pfilter"),
    function(object) {
        object@loglik
    }
)
