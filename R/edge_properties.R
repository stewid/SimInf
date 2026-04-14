## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
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

##' Convert an edge list with properties to a matrix
##'
##' A utility function to convert a data frame of edge properties
##' (e.g., movement rates, counts) into the specific matrix format
##' required for \code{ldata} in spatial models. This format allows
##' the C code to efficiently iterate over neighbors and their
##' associated properties.
##'
##' The function converts the \code{edges} data frame into a numeric matrix where:
##' \itemize{
##'   \item Each \strong{column} corresponds to a \code{to} node.
##'   \item Each column contains a sequence of blocks, one for each
##'   incoming edge.
##'   \item Each block starts with the \strong{zero-based index} of
##'   the \code{from} node.
##'   \item The subsequent rows in the block contain the property
##'   values (e.g., rate, count).
##'   \item Each block ends with a \strong{stop marker} (\code{-1}) in
##'   the first column.
##'   \item Unused cells in the matrix are filled with \code{NaN}.
##' }
##'
##' The \code{edges} data frame must contain columns \code{from} and
##' \code{to} with valid node indices (1-based). All edges pointing to
##' the same \code{to} node must be unique (no duplicate \code{from ->
##' to} pairs).
##'
##' @param edges A \code{data.frame} with properties assigned for each
##'     edge (\code{from} --> \code{to}). Must contain columns
##'     \code{from} and \code{to} (integer indices, 1-based) and any
##'     number of numeric property columns.
##' @param n_nodes The total number of nodes in the model. The
##'     resulting matrix will have \code{n_nodes} columns.
##' @return A numeric matrix with dimensions determined by the maximum
##'     number of incoming edges for any node. Columns correspond to
##'     \code{to} nodes.  Entries are \code{NaN} where no data exists.
##' @export
##' @examples
##' ## Define edge properties: from, to, rate, and count.
##' edges <- data.frame(
##'   from  = c(2, 3, 4, 1, 4, 5, 1, 3, 1, 3),
##'   to    = c(1, 1, 1, 2, 3, 3, 4, 4, 5, 5),
##'   rate  = c(0.2, 0.01, 0.79, 1, 0.2, 0.05, 0.2, 0.8, 0.2, 0.8),
##'   count = c(5, 5, 5, 50, 10, 10, 5, 5, 5, 5)
##' )
##'
##' ## Convert to matrix for 6 nodes.
##' ## Note: NaN values are printed as '.' for clarity.
##' mat <- edge_properties_to_matrix(edges, 6)
##' print(mat, na.print = ".")
##'
##' ## Interpretation of Column 1 (to node 1):
##' ## Row 1: Index of 'from' node (2-1 = 1), then stop marker (-1) at row 4.
##' ## Block 1: from=2 (idx 1), rate=0.2, count=5.
##' ## Block 2: from=3 (idx 2), rate=0.01, count=5.
##' ## Block 3: from=4 (idx 3), rate=0.79, count=5.
edge_properties_to_matrix <- function(edges, n_nodes) {
    ## Check that all values are numeric and finite.
    if (!all(sapply(edges, is.numeric)) ||
        !all(sapply(edges, is.finite))) {
        stop("Values in 'edges' must be numeric and finite.",
             call. = FALSE)
    }

    ## Check that all indices 'from -> to' are integers.
    if (any(!all(is_wholenumber(edges$from)),
            !all(is_wholenumber(edges$to)))) {
        stop("'edges' contain invalid 'from -> to' indices.",
             call. = FALSE)
    }

    edges$from <- as.integer(edges$from)
    edges$to <- as.integer(edges$to)

    ## Check that all indices 'from -> to' are valid.
    if (any(any(edges$from < 1L),
            any(edges$from > n_nodes),
            any(edges$to < 1L),
            any(edges$to > n_nodes),
            any(edges$from == edges$to))) {
        stop("'edges' contain invalid 'from -> to' indices.",
             call. = FALSE)
    }

    ## Check that all event pairs 'from -> to' are unique.
    if (any(tapply(edges$to, edges$from, anyDuplicated))) {
        stop("'edges' contain duplicated properties.",
             call. = FALSE)
    }

    ## Create a matrix with sequences starting with -1, indicating the
    ## stop condition for each item. Set all other values to
    ## 'NA'. Ensure that there is at least one stop value '-1' in each
    ## column after the values have been added to the matrix.
    n_to <- table(edges$to)
    m <- matrix(NaN,
                nrow = max(n_to) * (ncol(edges) - 1L) + 1L,
                ncol = n_nodes)

    ## Set stop values.
    m[1, ] <- -1
    m[cbind(n_to * (ncol(edges) - 1L) + 1L, as.integer(names(n_to)))] <- -1

    ## Determine the row index in the matrix to the first item
    ## ('from') of each sequence.
    i <- stats::ave(rep(1L, nrow(edges)), edges$to, FUN = cumsum)
    i <- (i - 1L) * (ncol(edges) - 1L) + 1L

    ## Set 'from', the first item in each sequence, to the zero-based
    ## index to 'from', i.e., decrease 'from' by 1.
    m[cbind(i, edges$to)] <- edges$from - 1

    ## Set all other values (except 'from' and 'to') in each sequence.
    j <- seq_len(ncol(edges))
    j <- j[!(colnames(edges)[j] %in% c("from", "to"))]
    for (k in seq_along(j)) {
        m[cbind(i + k, edges$to)] <- edges[, j[k], drop = TRUE]
    }

    m
}
