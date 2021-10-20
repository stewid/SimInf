## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
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

##' Box plot of number of individuals in each compartment
##'
##' Produce box-and-whisker plot(s) of the number of individuals in
##' each model compartment.
##' @param x The \code{model} to plot
##' @param compartments specify the names of the compartments to
##'     extract data from. The compartments can be specified as a
##'     character vector e.g. \code{compartments = c('S', 'I', 'R')},
##'     or as a formula e.g. \code{compartments = ~S+I+R} (see
##'     \sQuote{Examples}). Default (\code{compartments=NULL})
##'     includes all compartments.
##' @param index indices specifying the nodes to include when plotting
##'     data. Default \code{index = NULL} include all nodes in the
##'     model.
##' @param ... Additional arguments affecting the plot produced.
##' @aliases boxplot,SimInf_model-method
##' @export
##' @include SimInf_model.R
##' @importFrom graphics boxplot
##' @examples
##' ## Create an 'SIR' model with 10 nodes and initialise
##' ## it with 99 susceptible individuals and one infected
##' ## individual. Let the model run over 100 days.
##' model <- SIR(u0 = data.frame(S = rep(99, 10),
##'                              I = rep(1, 10),
##'                              R = rep(0, 10)),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the model and save the result.
##' result <- run(model)
##'
##' ## Create a boxplot that includes all compartments in all nodes.
##' boxplot(result)
##'
##' ## Create a boxplot that includes the S and I compartments in
##' ## nodes 1 and 2.
##' boxplot(result, ~S+I, 1:2)
setMethod(
    "boxplot",
    signature(x = "SimInf_model"),
    function(x, compartments = NULL, index = NULL, ...) {
        ## Remove the first two columns node and time
        boxplot(trajectory(x, compartments, index)[c(-1, -2)], ...)
    }
)

##' Scatterplot of number of individuals in each compartment
##'
##' A matrix of scatterplots with the number of individuals in each
##' compartment is produced. The \code{ij}th scatterplot contains
##' \code{x[,i]} plotted against \code{x[,j]}.
##' @param x The \code{model} to plot
##' @param compartments specify the names of the compartments to
##'     extract data from. The compartments can be specified as a
##'     character vector e.g. \code{compartments = c('S', 'I', 'R')},
##'     or as a formula e.g. \code{compartments = ~S+I+R} (see
##'     \sQuote{Examples}). Default (\code{compartments=NULL})
##'     includes all compartments.
##' @param index indices specifying the nodes to include when plotting
##'     data. Default \code{index = NULL} include all nodes in the
##'     model.
##' @param ... Additional arguments affecting the plot produced.
##' @export
##' @include SimInf_model.R
##' @importFrom graphics pairs
##' @examples
##' ## Create an 'SIR' model with 10 nodes and initialise
##' ## it with 99 susceptible individuals and one infected
##' ## individual. Let the model run over 100 days.
##' model <- SIR(u0 = data.frame(S = rep(99, 10),
##'                              I = rep(1, 10),
##'                              R = rep(0, 10)),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the model and save the result.
##' result <- run(model)
##'
##' ## Create a scatter plot that includes all compartments in all
##' ## nodes.
##' pairs(result)
##'
##' ## Create a scatter plot that includes the S and I compartments in
##' ## nodes 1 and 2.
##' pairs(result, ~S+I, 1:2)
setMethod(
    "pairs",
    signature(x = "SimInf_model"),
    function(x, compartments = NULL, index = NULL, ...) {
        ## Remove the first two columns node and time
        pairs(trajectory(x, compartments, index)[c(-1, -2)], ...)
    }
)

init_plot_node_index <- function(model, index) {
    index <- check_node_index_argument(model, index)
    if (is.null(index))
        index <- seq_len(n_nodes(model))
    index
}

init_plot_line_type <- function(lty, compartments, each) {
    if (is.null(compartments)) {
        n <- 1
    } else {
        n <- length(compartments)
    }

    if (is.null(lty)) {
        lty <- seq_len(n)
    } else {
        lty <- rep(lty, length.out = n)
    }
    rep(lty, each = each)
}

init_plot_color <- function(col, compartments, each) {
    if (is.null(compartments)) {
        n <- 1
    } else {
        n <- length(compartments)
    }

    if (is.null(col)) {
        if (n > 9) {
            col <- rainbow(n)
        } else if (n > 1) {
            col <- rep(c("#e41a1c", "#377eb8", "#4daf4a",
                         "#984ea3", "#ff7f00", "#ffff33",
                         "#a65628", "#f781bf", "#999999"),
                       length.out = n)
        } else {
            col <- "black"
        }
    } else {
        col <- rep(col, length.out = n)
    }
    rep(col, each = each)
}

init_plot_range <- function(range) {
    if (identical(range, FALSE))
        return(range)

    if (any(!is.numeric(range),
            !identical(length(range), 1L),
            range < 0, range > 1)) {
        stop("'range' must be FALSE or a value between 0 and 1.",
             call. = FALSE)
    }

    (1 - range) / 2
}

init_plot_prevalence_data <- function(model, compartments,
                                      level, index, range) {
    index <- init_plot_node_index(model, index)
    range <- init_plot_range(range)

    ## Create a matrix with one row for each line in the plot.
    if (level < 3) {
        y <- prevalence(model, compartments, level, index, "matrix")
        lower <- NULL
        upper <- NULL
        each <- 1
    } else if (identical(range, FALSE) || identical(length(index), 1L)) {
        y <- prevalence(model, compartments, level, index, "matrix")
        lower <- NULL
        upper <- NULL
        each <- length(index)
    } else {
        y <- apply(
            prevalence(model, compartments, level, index, "matrix"),
            2, quantile, probs = c(range, 0.5, 1 - range))

        ## Matrices for quantile ranges and median.
        j <- seq_len(ncol(y))
        lower <- y[1, j, drop = FALSE]
        upper <- y[3, j, drop = FALSE]
        y <- y[2, j, drop = FALSE]
        each <- 1
    }

    list(lower        = lower,
         y            = y,
         upper        = upper,
         each         = each,
         compartments = NULL)
}

init_plot_trajectory_data <- function(model, compartments, index, range) {
    index <- init_plot_node_index(model, index)
    range <- init_plot_range(range)

    ## Create a matrix with one row for each line in the plot.
    y <- list()
    for (j in seq_len(length(compartments$rhs))) {
        for (compartment in names(compartments$rhs[[j]])) {
            if (identical(range, FALSE) || identical(length(index), 1L)) {
                y[[length(y) + 1]] <-
                    trajectory(model, compartment, index, "matrix")
            } else {
                y[[length(y) + 1]] <- apply(
                    trajectory(model, compartment, index, "matrix"),
                    2, quantile, probs = c(range, 0.5, 1 - range))
            }

            names(y)[length(y)] <- compartment
        }
    }

    compartments <- names(y)

    if (identical(range, FALSE) || identical(length(index), 1L)) {
        lower <- NULL
        upper <- NULL
        ## Combine matrices for each comparment.
        y <- do.call("rbind", y)
        each <- length(index)
    } else {
        ## Matrices for quantile ranges and median.
        lower <- do.call("rbind", lapply(y, function(x) x[1, ]))
        upper <- do.call("rbind", lapply(y, function(x) x[3, ]))
        y <- do.call("rbind", lapply(y, function(x) x[2, ]))
        each <- 1
    }

    list(lower        = lower,
         y            = y,
         upper        = upper,
         each         = each,
         compartments = compartments)
}

##' Determine if the 'compartments' expression contains a lhs.
##' @noRd
compartments_has_lhs <- function(compartments) {
    if (is(compartments, "formula")) {
        compartments <- as.character(compartments)
        if (identical(length(compartments), 3L))
            return(TRUE)
    }

    FALSE
}

init_plot_argv <- function(model, compartments, pd, type, lwd, ...) {
    argv <- list(...)
    argv$type <- type
    argv$lwd <- lwd

    if (is.null(argv$ylab)) {
        if (isTRUE(compartments_has_lhs(compartments))) {
            argv$ylab <- deparse(compartments)
        } else {
            argv$ylab <- "Value"
        }
    }

    ## Settings for the y-axis.
    if (is.null(argv$ylim)) {
        if (is.null(pd$upper)) {
            argv$ylim <- c(0, max(pd$y))
        } else {
            argv$ylim <- c(0, max(pd$upper))
        }
    }

    ## Settings for the x-axis
    if (is.null(names(model@tspan))) {
        argv$x <- model@tspan
        if (is.null(argv$xlab))
            argv$xlab <- "Time"
    } else {
        argv$x <- as.Date(names(model@tspan))
        if (is.null(argv$xlab))
            argv$xlab <- "Date"
    }

    argv
}

plot_data <- function(pd, argv, lty, col, frame.plot, legend) {
    ## Plot lines
    for (i in seq_len(dim(pd$y)[1])) {
        argv$y <- pd$y[i, ]
        argv$col <- col[i]
        argv$lty <- lty[i]

        if (i == 1) {
            argv$frame.plot <- frame.plot
            do.call(plot, argv)
            argv$frame.plot <- NULL
        } else {
            do.call(lines, argv)
        }

        if (!is.null(pd$lower) && !is.null(pd$upper)) {
            if (argv$type == "s") {
                x <- c(rep(argv$x, each = 2)[-1],
                       rep(rev(argv$x), each = 2)[-1])
                y <- c(rep(pd$upper[i, ], each = 2)[-2 * ncol(pd$upper)],
                       rep(rev(pd$lower[i, ]), each = 2)[-2 * ncol(pd$lower)])
            } else {
                x <- c(argv$x, rev(argv$x))
                y <- c(pd$upper[i, ], rev(pd$lower[i, ]))
            }

            polygon(x = x, y = y, border = NA,
                    col = adjustcolor(col[i], alpha.f = 0.1))
        }
    }

    ## Add the legend below plot. The default legend is the names of
    ## the compartments.
    if (isTRUE(legend) && !is.null(pd$compartments)) {
        ## Determine the size of the legend.
        lgd <- legend("top", lty = unique(lty), col = unique(col),
                      bty = "n", horiz = TRUE, legend = pd$compartments,
                      lwd = argv$lwd, xpd = TRUE, plot = FALSE)

        ## Determine the y-position to place the legend one line above
        ## the plot.
        y <- par("cin")[2] * par("cex") * par("lheight")
        y <- diff(grconvertY(c(0, y), "inches", "npc"))
        y <- grconvertY(1 + y, "npc", "user")

        ## Add the legend to the plot.
        legend(lgd$rect$left, y, lty = unique(lty), col = unique(col),
               bty = "n", horiz = TRUE, legend = pd$compartments,
               lwd = argv$lwd, xpd = TRUE)
    }
}

##' @importFrom graphics contour
##' @importFrom graphics lines
##' @importFrom graphics rug
##' @importFrom MASS kde2d
##' @importFrom stats density
##' @noRd
plot_density <- function(x, ...) {
    if (ncol(x) > 1) {
        pairs(x,
              diag.panel = function(x, ...) {
                  usr <- par("usr")
                  on.exit(par(usr))
                  par(usr = c(usr[1:2], 0, 1.5))
                  d <- density(x, bw = "SJ-ste")
                  d$y <- d$y / max(d$y)
                  lines(d, ...)
                  rug(x)
              },
              lower.panel = function(x, y, ...) {
                  d <- kde2d(x, y)
                  contour(d, add = TRUE, drawlabels = FALSE, ...)
              }, ...)
    } else {
        plot(density(x, bw = "SJ-ste"), main = "", xlab = colnames(x), ...)
        rug(x)
    }
}

plot_trace <- function(x, i, j, ...) {
    mfrow <- switch(min(length(j), 10),
                    c(1, 1),
                    c(1, 2),
                    c(2, 2),
                    c(2, 2),
                    c(3, 2),
                    c(3, 2),
                    c(3, 3),
                    c(3, 3),
                    c(3, 3),
                    stop("To many variables to plot."))
    opar <- par(mfrow = mfrow)
    on.exit(par(opar), add = TRUE)

    for (k in j) {
        plot(x[i, k] ~ i,
             xlab = "Iterations",
             ylab = colnames(x)[k],
             type = "l",
             main = sprintf("Trace of %s", colnames(x)[k]),
             ...)
    }
}

##' Display the outcome from a simulated trajectory
##'
##' Plot either the median and the quantile range of the counts in all
##' nodes, or plot the counts in specified nodes.
##' @param x The \code{model} to plot.
##' @template plot-y-param
##' @template plot-level-param
##' @template plot-index-param
##' @template plot-range-param
##' @template plot-type-param
##' @template plot-lwd-param
##' @template plot-frame-param
##' @template plot-legend-param
##' @param ... Other graphical parameters that are passed on to the
##'     plot function.
##' @rdname plot
##' @aliases plot,SimInf_model-method
##' @export
##' @include SimInf_model.R
##' @importFrom graphics grconvertY
##' @importFrom graphics legend
##' @importFrom graphics lines
##' @importFrom graphics par
##' @importFrom graphics plot
##' @importFrom graphics polygon
##' @importFrom grDevices adjustcolor
##' @importFrom grDevices rainbow
##' @examples
##' ## Create an 'SIR' model with 100 nodes and initialise
##' ## it with 990 susceptible individuals and 10 infected
##' ## individuals in each node. Run the model over 100 days.
##' model <- SIR(u0 = data.frame(S = rep(990, 100),
##'                              I = rep(10, 100),
##'                              R = rep(0, 100)),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the model and save the result.
##' result <- run(model)
##'
##' ## Plot the median and interquartile range of the number
##' ## of susceptible, infected and recovered individuals.
##' plot(result)
##'
##' ## Plot the median and the middle 95\% quantile range of the
##' ## number of susceptible, infected and recovered individuals.
##' plot(result, range = 0.95)
##'
##' ## Plot the median and interquartile range of the  number
##' ## of infected individuals.
##' plot(result, "I")
##'
##' ## Use the formula notation instead to plot the median and
##' ## interquartile range of the number of infected individuals.
##' plot(result, ~I)
##'
##' ## Plot the number of susceptible, infected
##' ## and recovered individuals in the first
##' ## three nodes.
##' plot(result, index = 1:3, range = FALSE)
##'
##' ## Use plot type line instead.
##' plot(result, index = 1:3, range = FALSE, type = "l")
##'
##' ## Plot the number of infected individuals in the first node.
##' plot(result, "I", index = 1, range = FALSE)
##'
##' ## Plot the proportion of infected individuals (cases)
##' ## in the population.
##' plot(result, I ~ S + I + R)
##'
##' ## Plot the proportion of nodes with infected individuals.
##' plot(result, I ~ S + I + R, level = 2)
##'
##' ## Plot the median and interquartile range of the proportion
##' ## of infected individuals in each node
##' plot(result, I ~ S + I + R, level = 3)
##'
##' ## Plot the proportion of infected individuals in the first
##' ## three nodes.
##' plot(result, I ~ S + I + R, level = 3, index = 1:3, range = FALSE)
setMethod(
    "plot",
    signature(x = "SimInf_model", y = "ANY"),
    function(x, y, level = 1, index = NULL, range = 0.5, type = "s",
             lwd = 2, frame.plot = FALSE, legend = TRUE, ...) {
        if (missing(y))
            y <- NULL

        if (isTRUE(compartments_has_lhs(y))) {
            pd <- init_plot_prevalence_data(x, y, level, index, range)
        } else {
            compartments <- match_compartments(compartments = y,
                                               ok_combine = TRUE,
                                               ok_lhs = FALSE,
                                               U = rownames(x@S),
                                               V = rownames(x@v0))

            pd <- init_plot_trajectory_data(x, compartments, index, range)
        }

        argv <- init_plot_argv(x, y, pd, type, lwd, ...)
        lty <- init_plot_line_type(argv$lty, pd$compartments, pd$each)
        col <- init_plot_color(argv$col, pd$compartments, pd$each)

        plot_data(pd, argv, lty, col, frame.plot, legend)

        invisible(NULL)
    }
)

##' Display the ABC posterior distribution
##'
##' @param x The \code{SimInf_abc} object to plot.
##' @param y The generation to plot. The default is to display the
##'     last generation.
##' @param ... Additional arguments affecting the plot.
##' @aliases plot,SimInf_abc-method
##' @export
##' @include abc.R
setMethod(
    "plot",
    signature(x = "SimInf_abc"),
    function(x, y, ...) {
        if (missing(y))
            y <- length(x@x)
        y <- as.integer(y)
        if (length(y) != 1) {
            stop("Can only select one generation to plot.",
                 call. = FALSE)
        }

        plot_density(t(x@x[[y]]), ...)

        invisible(NULL)
    }
)

##' Display the PMCMC posterior distribution
##'
##' Display the (approximate) posterior distributions obtained from
#'      fitting a particle Markov chain Monte Carlo algorithm, or the
#'      corresponding trace plots.
##' @param x The \code{SimInf_pmcmc} object to plot.
##' @param y The trace of all variables and logPost are plotted when
##'     \code{y = "trace"} or \code{y = ~trace}, else the posterior
##'     distributions are plotted. Default is to plot the posterier
##'     distributions.
##' @param start The start iteration to remove some burn-in
##'     iterations. Default is \code{start = 1}.
##' @param thin keep every \code{thin} iteration after the
##'     \code{start} iteration. Default is \code{thin = 1}, i.e., keep
##'     every iteration.
##' @param ... Additional arguments affecting the plot.
##' @aliases plot,SimInf_pmcmc-method
##' @export
##' @include pmcmc.R
setMethod(
    "plot",
    signature(x = "SimInf_pmcmc"),
    function(x, y, start = 1, thin = 1, ...) {
        if (missing(y))
            y <- "density"
        if (is(y, "formula") && identical(as.character(y), c("~", "trace")))
            y <- "trace"
        y <- as.character(y)

        start <- as.integer(start)
        if (any(length(start) != 1, any(start < 1)))
            stop("'start' must be an integer >= 1.", call. = FALSE)

        thin <- as.integer(thin)
        if (any(length(thin) != 1, any(thin < 1)))
            stop("'thin' must be an integer >= 1.", call. = FALSE)

        i <- seq(from = start, to = length(x), by = thin)
        j <- seq(from = 5, by = 1, length.out = length(x@pars))

        if (identical(y, "trace")) {
            ## Include the logPost column.
            plot_trace(x@chain, i, c(j, 1), ...)
        } else {
            plot_density(x@chain[i, j, drop = FALSE], ...)
        }

        invisible(NULL)
    }
)
