## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
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

##' Create a distance matrix between nodes for spatial models
##'
##' Calculate the Euclidean distances between all pairs of nodes based
##' on their projected coordinates (\code{x}, \code{y}). Distances
##' greater than the specified \code{cutoff} are excluded from the
##' result (stored as zeros in the sparse matrix).
##'
##' The result is a symmetric sparse matrix (\code{dgCMatrix}) where
##' the element \code{d[i, j]} contains the distance between node
##' \code{i} and node \code{j} if it is less than or equal to
##' \code{cutoff}, and \code{0} otherwise.
##'
##' @param x Numeric vector of projected x coordinates for each node.
##' @param y Numeric vector of projected y coordinates for each node.
##' @param cutoff Numeric scalar. The maximum distance to include in
##'     the matrix. Pairs of nodes farther apart than this value are
##'     excluded from the sparse structure (stored as zeros).
##' @param min_dist Numeric scalar. The value to use for the distance
##'     between two nodes if their coordinates are identical (distance
##'     = 0).  This prevents division by zero errors in downstream
##'     calculations (e.g., inverse distance weighting). If
##'     \code{NULL} (default) and identical coordinates are found, an
##'     error is raised.
##' @param na_fail Logical. If \code{TRUE} (default), missing values
##'     (\code{NA}) in \code{x} or \code{y} will raise an error. If
##'     \code{FALSE}, distances involving missing coordinates are set
##'     to zero.
##' @return A symmetric sparse matrix of class
##'     \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}.  Non-zero
##'     entries represent distances \eqn{\le} \code{cutoff}; entries
##'     outside the cutoff are implicitly zero.
##' @export
##' @examples
##' ## Generate a 10 x 10 grid of nodes separated by 100m.
##' nodes <- expand.grid(
##'   x = seq(from = 0, to = 900, by = 100),
##'   y = seq(from = 0, to = 900, by = 100)
##' )
##' plot(nodes, main = "Node Grid", asp = 1)
##'
##' ## Calculate distances with a 300m cutoff.
##' ## Only neighbors within 300m will have non-zero entries.
##' d <- distance_matrix(
##'   x = nodes$x,
##'   y = nodes$y,
##'   cutoff = 300
##' )
##'
##' ## Inspect the sparse matrix structure.
##' ## Note: The matrix is symmetric and the diagonal is zero.
##' d[1:10, 1:10]
##'
##' ## Count the number of neighbors for the first node.
##' sum(d[1, ] > 0)
distance_matrix <- function(x, y, cutoff, min_dist = NULL, na_fail = TRUE) {
    .Call(SimInf_distance_matrix,
          as.numeric(x),
          as.numeric(y),
          as.numeric(cutoff),
          as.numeric(min_dist),
          isTRUE(na_fail))
}

##' Add spatial coupling information to local data
##'
##' A utility function to augment local model parameters
##' (\code{ldata}) with spatial coupling data from neighboring nodes.
##'
##' The function calculates distances between nodes based on projected
##' coordinates (\code{x}, \code{y}) and appends neighbor data to the
##' \code{ldata} matrix for each node.
##'
##' \strong{Output Format:}
##' The returned matrix has the same number of rows as the input
##' \code{ldata}.  The columns are organized as follows:
##' \itemize{
##'   \item \strong{Local Parameters}: The first \eqn{n} columns
##'     correspond to the original local parameters passed in
##'     \code{ldata}.
##'   \item \strong{Neighbor Pairs}: Following the local parameters,
##'     the data is stored as pairs of columns: \code{(neighbor_index,
##'     distance)}.  The \code{neighbor_index} is a zero-based index
##'     of the neighbor node.  The \code{distance} is the Euclidean
##'     distance to that neighbor.
##'   \item \strong{Stop Marker}: Each node's neighbor list is
##'     terminated by a pair \code{(-1, 0)} in the position where the
##'     next \code{neighbor_index} would appear. This marker appears
##'     exactly once per node, immediately after the last neighbor.
##' }
##'
##' @param x Numeric vector of projected x coordinates for each node.
##' @param y Numeric vector of projected y coordinates for each node.
##' @param cutoff Numeric scalar. The maximum distance for considering
##'     two nodes as neighbors. Pairs of nodes farther apart than this
##'     value are excluded from the neighbor data.
##' @template ldata-param
##' @param min_dist Numeric scalar. The value to use for the distance
##'     between two nodes if their coordinates are identical (distance
##'     = 0).  This prevents division by zero errors in downstream
##'     calculations (e.g., inverse distance weighting). If
##'     \code{NULL} (default) and identical coordinates are found, an
##'     error is raised.
##' @param na_fail Logical. If \code{TRUE} (default), missing values
##'     (\code{NA}) in \code{x} or \code{y} will raise an error. If
##'     \code{FALSE}, distances involving missing coordinates are set
##'     to zero.
##' @return A numeric matrix with the same number of rows as
##'     \code{ldata}, but with additional columns containing the
##'     neighbor indices and distances as described in the
##'     \sQuote{Output Format} section.
##' @seealso \code{\link{distance_matrix}} for computing pairwise
##'     distances between nodes, and
##'     \code{\link{edge_properties_to_matrix}} for a similar utility
##'     that converts edge properties to a matrix format.
##' @export
##' @examples
##' ## Generate a 5 x 5 grid of nodes separated by 1000m.
##' nodes <- expand.grid(
##'   x = seq(from = 0, to = 4000, by = 1000),
##'   y = seq(from = 0, to = 4000, by = 1000)
##' )
##'
##' ## Create local data with one parameter per node.
##' ldata <- matrix(0.1, nrow = 1, ncol = nrow(nodes))
##'
##' ## Add spatial coupling with a 2500m cutoff.
##' ldata_augmented <- add_spatial_coupling_to_ldata(
##'   x = nodes$x,
##'   y = nodes$y,
##'   cutoff = 2500,
##'   ldata = ldata
##' )
##'
##' ## Inspect the result for the first node.
##' ldata_augmented[, 1]
add_spatial_coupling_to_ldata <- function(x, y, cutoff, ldata = NULL,
                                          min_dist = NULL, na_fail = TRUE) {
    distance <- distance_matrix(x = x,
                                y = y,
                                cutoff = cutoff,
                                min_dist = min_dist,
                                na_fail = na_fail)

    if (is.null(ldata)) {
        ldata <- matrix(data = numeric(), nrow = 0, ncol = ncol(distance))
    } else {
        ldata <- init_data_matrix(ldata)

        if (ncol(ldata) != ncol(distance)) {
            stop("Number of nodes in 'ldata' and coordinates must match.",
                 call. = FALSE)
        }
    }

    .Call(SimInf_ldata_sp, ldata, distance, 1L)
}
