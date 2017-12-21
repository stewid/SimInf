## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2017  Stefan Engblom
## Copyright (C) 2015 - 2017  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' Create a distance matrix between nodes for spatial models
##'
##' Calculate the euclidian distances beween coordinates for all
##' coordinates within the cutoff.
##' @param x Projected x coordinate
##' @param y Projected y coordinate
##' @param cutoff The distance cutoff
##' @param min_dist The minimum distance to separate two nodes.  If
##'     the coordinates for two nodes are identical, the min_dist must
##'     be assigned or an error is raised.  Default is NULL i.e. to
##'     raise an error.
##' @return \code{dgCMatrix}
##' @export
##' @importFrom methods new
##' @examples
##' ## Generate a grid 10 x 10 and place one node in each cell
##' ## separated by 100m.
##' nodes <- expand.grid(x = (0:9) * 100, y = (0:9) * 100)
##' plot(y ~ x, nodes)
##'
##' ## Define the cutoff to only include neighbors within 300m.
##' d <- distance_matrix(x = nodes$x, y = nodes$y, cutoff = 301)
##'
##' ## View the first 10 rows and columns in the distance matrix
##' d[1:10, 1:10]
distance_matrix <- function(x, y, cutoff, min_dist = NULL)
{
    if (!is.null(min_dist)) {
        if (any(!is.numeric(min_dist),
                !identical(length(min_dist), 1L),
                min_dist[1] <= 0))
        {
            stop("Invalid 'min_dist' argument. Please provide 'min_dist' > 0.")
        }
    }

    m <- lapply(seq_len(length(x)), function(i) {
        x0 <- x
        y0 <- y
        x1 <- x0[i]
        y1 <- y0[i]

        ## Calculate euclidian distance
        d <- sqrt((x0 - x1)^2 + (y0 - y1)^2)

        ## Determine which indices are closer than cutoff
        row_ind <- which(d < cutoff)

        ## Drop current i
        row_ind <- row_ind[row_ind != i]

        d <- d[row_ind]

        if (any(d == 0)) {
            if (is.null(min_dist))
                stop("Identical coordinates. Please provide a minimum distance.")
            d <- sapply(d, max, min_dist)
        }

        if (!is.null(min_dist))
            d <- sapply(d, max, min_dist)

        ## Make row indices 0-based
        list(row_ind = row_ind - 1, d = d)
    })

    ## Create vectors for all distances, row indices and column indices.
    d <- as.numeric(unlist(sapply(m, "[[", "d")))
    row_ind <- as.integer(unlist(sapply(m, "[[", "row_ind")))
    col_ind <- as.integer(c(0, cumsum(sapply(m, function(x) length(x$row_ind)))))

    ## Create a new sparse matrix
    new("dgCMatrix", x = d, i = row_ind, p = col_ind, Dim = rep(length(x), 2))
}
