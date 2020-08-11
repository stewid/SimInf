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
##' ## Create a boxplot
##' boxplot(result)
setMethod(
    "boxplot",
    signature(x = "SimInf_model"),
    function(x, ...) {
        ## Remove the first two columns node and time
        boxplot(trajectory(x)[c(-1, -2)], ...)
    }
)

##' Scatterplot of number of individuals in each compartment
##'
##' A matrix of scatterplots with the number of individuals in each
##' compartment is produced. The \code{ij}th scatterplot contains
##' \code{x[,i]} plotted against \code{x[,j]}.
##' @param x The \code{model} to plot
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
##' ## Create a scatter plot
##' pairs(result)
setMethod(
    "pairs",
    signature(x = "SimInf_model"),
    function(x, ...) {
        ## Remove the first two columns node and time
        pairs(trajectory(x)[c(-1, -2)], ...)
    }
)

init_plot_compartments <- function(x, compartments) {
    ## Determine the compartments to include in the plot
    if (is.null(compartments))
        compartments <- rownames(x@S)
    if (!(all(compartments %in% rownames(x@S))))
        stop("'compartments' must exist in the model.", call. = FALSE)
    match(compartments, rownames(x@S))
}

init_plot_node <- function(x, node) {
    node <- check_node_index_argument(x, node)
    if (is.null(node))
        node <- seq_len(n_nodes(x))
    node
}

init_plot_line_type <- function(lty, compartments, m) {
    if (is.null(lty)) {
        lty <- seq_len(length(compartments))
    } else {
        lty <- rep(lty, length.out = length(compartments))
    }
    rep(lty, length.out = dim(m)[1])
}

init_plot_color <- function(col, compartments, m) {
    if (is.null(col)) {
        if (length(compartments) > 9) {
            col <- rainbow(length(compartments))
        } else if (length(compartments) > 1) {
            col <- rep(c("#e41a1c", "#377eb8", "#4daf4a",
                         "#984ea3", "#ff7f00", "#ffff33",
                         "#a65628", "#f781bf", "#999999"),
                       length.out = length(compartments))
        } else {
            col <- "black"
        }
    } else {
        col <- rep(col, length.out = length(compartments))
    }
    rep(col, length.out = dim(m)[1])
}

init_plot_type <- function(type) {
    if (is.null(type))
        type <- "l"
    type
}

init_plot_line_width <- function(lwd) {
    if (is.null(lwd))
        lwd <- 2
    lwd
}

##' Display the outcome from a simulated trajectory
##'
##' Plot either the median and the quantile range of the counts in all
##' nodes, or plot the counts in specified nodes.
##' @param x The \code{model} to plot.
##' @param compartments Character vector with the compartments in the
##'     model to include in the plot. Default is \code{NULL}
##'     i.e. include all compartments in the model.
##' @param node Indices specifying the nodes to include when plotting
##'     data. Plot one line for each node. Default (\code{node =
##'     NULL}) is to extract data from all nodes and plot the median
##'     count for the specified compartments.
##' @param range Show the quantile range of the count in each
##'     compartment. Default is to show the interquartile range
##'     i.e. the middle 50\% of the count in transparent color. The
##'     median value is shown in the same color. Use \code{range =
##'     0.95} to show the middle 95\% of the count. To display
##'     individual lines for each node, specify \code{range = FALSE}.
##' @param ... Additional arguments affecting the plot produced.
##' @rdname plot
##' @aliases plot,SimInf_model-method
##' @export
##' @include SimInf_model.R
##' @importFrom graphics legend
##' @importFrom graphics lines
##' @importFrom graphics par
##' @importFrom graphics plot
##' @importFrom graphics polygon
##' @importFrom graphics title
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
##' plot(result, compartments = "I")
##'
##' ## Plot the number of susceptible, infected
##' ## and recovered individuals in the first
##' ## three nodes.
##' plot(result, node = 1:3, range = FALSE)
##'
##' ## Use a stair step plot type instead.
##' plot(result, node = 1:3, range = FALSE, type = "s")
##'
##' ## Plot the number of infected individuals in the first node.
##' plot(result, compartments = "I", node = 1, range = FALSE)
setMethod(
    "plot",
    signature(x = "SimInf_model"),
    function(x, compartments = NULL, node = NULL, range = 0.5, ...) {
        if (identical(dim(x@U), c(0L, 0L))) {
            stop("Please run the model first, the 'U' matrix is empty.",
                 call. = FALSE)
        }

        argv <- list(...)

        compartments <- init_plot_compartments(x, compartments)
        node <- init_plot_node(x, node)

        ## Create a matrix with one row for each line in the plot.
        if (identical(range, FALSE)) {
            ## Extract subset of data from U
            i <- rep(compartments, length(node))
            i <- i + rep((node - 1) * Nc(x), each = length(compartments))
            m <- x@U[i, seq_len(ncol(x@U)), drop = FALSE]
        } else {
            ## Check range argument
            if (any(!is.numeric(range), !identical(length(range), 1L),
                    range < 0, range > 1)) {
                stop("'range' must be FALSE or a value between 0 and 1.",
                     call. = FALSE)
            }
            range <- (1 - range) / 2

            m <- matrix(0, nrow = length(compartments),
                        ncol = length(x@tspan))

            ## Matrices for quantile range
            mu <- m
            ml <- m

            for (i in seq_len(length(compartments))) {
                k <- seq(from = compartments[i], to = dim(x@U)[1],
                         by = Nc(x))
                u <- apply(x@U[k[node], seq_len(ncol(x@U)), drop = FALSE],
                           2,
                           quantile,
                           probs = c(range, 0.5, 1 - range))
                ml[i, ] <- u[1, ]
                m[i, ] <- u[2, ]
                mu[i, ] <- u[3, ]
            }

            range <- TRUE
        }

        lty <- init_plot_line_type(argv$lty, compartments, m)
        col <- init_plot_color(argv$col, compartments, m)
        argv$type <- init_plot_type(argv$type)
        argv$lwd <- init_plot_line_width(argv$lwd)

        ## Settings for the y-axis.
        argv$ylab <- "N"
        if (isTRUE(range)) {
            argv$ylim <- c(0, max(mu))
        } else {
            argv$ylim <- c(0, max(m))
        }

        ## Settings for the x-axis
        if (is.null(names(x@tspan))) {
            xx <- x@tspan
            argv$xlab <- "Time"
        } else {
            xx <- as.Date(names(x@tspan))
            argv$xlab <- "Date"
        }

        savepar <- par(mar = c(2, 4, 1, 1), oma = c(4, 1, 0, 0),
                       xpd = TRUE)
        on.exit(par(savepar))

        ## Plot lines
        for (i in seq_len(dim(m)[1])) {
            argv$x <- xx
            argv$y <- m[i, ]
            argv$col <- col[i]
            argv$lty <- lty[i]
            if (i == 1) {
                do.call(plot, argv)
                title(xlab = argv$xlab, outer = TRUE, line = 0)
            } else {
                do.call(lines, argv)
            }
            if (isTRUE(range)) {
                polygon(x = c(xx, rev(xx)), y = c(mu[i, ], rev(ml[i, ])),
                        col = adjustcolor(col[i], alpha.f = 0.1),
                        border = NA)
            }
        }

        ## Add the legend below plot. The default legend is the names
        ## of the compartments.
        if (is.null(argv$legend))
            argv$legend <- rownames(x@S)[compartments]
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
            mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        legend("bottom", inset = c(0, 0),
               lty = lty[seq_len(length(compartments))],
               col = col[seq_len(length(compartments))],
               bty = "n", horiz = TRUE, legend = argv$legend,
               lwd = argv$lwd)
    }
)

##' Display the ABC posterior distribution
##'
##' @param x The \code{SimInf_abc} object to plot.
##' @param y The generation to plot. The default is to display the
##'     last generation.
##' @param ... Additional arguments affecting the plot.
##' @aliases plot,SimInf_abc-method
##' @importFrom graphics contour
##' @importFrom graphics lines
##' @importFrom MASS bandwidth.nrd
##' @importFrom MASS kde2d
##' @importFrom stats density
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

        if (length(x@pars) > 1) {
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
        } else {
            plot(density(x@x[[y]], bw = "SJ-ste"), main = "",
                 xlab = rownames(x@x[[y]]), ...)
        }
    }
)
