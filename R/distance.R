## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2025 Stefan Widgren
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

##' Create a distance matrix between nodes for spatial models
##'
##' Calculate the euclidian distances beween coordinates for all
##' coordinates within the cutoff.
##' @param x Projected x coordinate
##' @param y Projected y coordinate
##' @param cutoff The distance cutoff
##' @param min_dist The minimum distance to separate two nodes.  If
##'     the coordinates for two nodes are identical, the min_dist must
##'     be assigned or an error is raised.  Default is \code{NULL},
##'     i.e., to raise an error.
##' @param na_fail A logical indicating whether missing values in
##'     \code{x} or \code{y} should raise an error or assign zero to
##'     all distances involving missing values.  Default is
##'     \code{TRUE}, i.e., to raise an error.
##' @return \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}
##' @export
##' @examples
##' ## Generate a grid 10 x 10 and place one node in each cell
##' ## separated by 100m.
##' nodes <- expand.grid(x = (0:9) * 100, y = (0:9) * 100)
##' plot(y ~ x, nodes)
##'
##' ## Define the cutoff to only include neighbors within 300m.
##' d <- distance_matrix(x = nodes$x, y = nodes$y, cutoff = 300)
##'
##' ## View the first 10 rows and columns in the distance matrix
##' d[1:10, 1:10]
distance_matrix <- function(x, y, cutoff, min_dist = NULL, na_fail = TRUE) {
    .Call(SimInf_distance_matrix,
          as.numeric(x),
          as.numeric(y),
          as.numeric(cutoff),
          as.numeric(min_dist),
          isTRUE(na_fail))
}

check_metric <- function(metric) {
    check_integer_arg(metric)
    metric <- as.integer(metric)
    if (any(length(metric) != 1L, any(metric < 0L), any(metric > 2L)))
        stop("'metric' must be an integer 0, 1, or 2.", call. = FALSE)
    metric
}

##' Add information about spatial coupling between nodes to 'ldata'
##'
##' A utility function to combine local model parameters
##' (\code{ldata}) and spatial coupling to other nodes and add the
##' result to \code{ldata}.
##'
##' Format for ldata: the first n indicies (1, 2, ..., n) are the
##' local model parameters, i.e, the indata to the function. They are
##' followed by the neighbor data, pairs of (index, value) and then a
##' stop pair (-1, 0) where 'index' is the zero-based index to the
##' neighbor and the value is determined by the metric argument.
##' @param ldata Matrix with local model parameters for each node.
##' @param distance Sparse matrix
##'     (\code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}) with
##'     distances between nodes, see \code{\link{distance_matrix}}.
##' @param metric The type of value to calculate for each neighbor: 0)
##'     degree, 1) distance, and 2) 1 / (distance * distance). Default
##'     is \code{1}.
##' @return matrix
##' @export
##' @examples
##' d <- distance_matrix(nodes$x, nodes$y, cutoff = 5000)
##' ldata <- rbind(nodes$x, nodes$y)
##' ldata <- ldata_add_spatial_coupling(ldata, d)
ldata_add_spatial_coupling <- function(ldata, distance, metric = 1) {
    metric <- check_metric(metric)
    .Call(SimInf_ldata_sp, ldata, distance, metric)
}
