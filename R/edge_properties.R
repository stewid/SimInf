## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2023 Stefan Widgren
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
##' A utility function to facilitate preparing edge properties for
##' \code{ldata} in a model.
##'
##' The edge properties will be converted to a matrix where each row
##' in \code{edges} will become a sequence of (index, value_1,
##' value_2, ..., value_n) where 'index' is the zero-based index of
##' the \code{from} node. The reason for a zero-based index is to
##' facilitate it's usage in C code. The sequence will be added to the
##' 'to' column in the matrix. There will always be at least one stop
##' value=-1 in each column. All other values in the matrix will be
##' set to \code{NaN}. See \sQuote{Examples}.
##' @param edges a \code{data.frame} with properties assigned for each
##'     edge 'from' --> 'to', for example, weight or count. The
##'     \code{data.frame} must contain the columns '\code{from}' and
##'     '\code{to}' with valid indices to the nodes (1 <= index <=
##'     n_nodes).
##' @param n_nodes the total number of nodes in the model. The
##'     resulting matrix will have the number of columns equal to
##'     \code{n_nodes}.
##' @return a numeric matrix with the number of rows equal to
##'     \code{max(table(edges$to)) * (ncol(edges) - 1) + 1} and the
##'     number of columns equal to \code{n_nodes}.
##' @export
##' @examples
##' ## Let us consider the following edge properties.
##' edges <- data.frame(
##'     from  = c(  2,    3,     4,  1,   4,    5,   1,   3,   1,   3),
##'     to    = c(  1,    1,     1,  2,   3,    3,   4,   4,   5,   5),
##'     rate  = c(0.2, 0.01,  0.79,  1, 0.2, 0.05, 0.2, 0.8, 0.2, 0.8),
##'     count = c(  5,    5,     5, 50,  10,   10,   5,   5,   5,   5))
##'
##' ## Converting the edge properties into a matrix
##' edge_properties_to_matrix(edges, 6)
##'
##' ## Gives the following output. The first column contains first the
##' ## properties for the edge from = 2 --> to = 1, where the first
##' ## row is the zero-based index of from, i.e., 1. The second row
##' ## contains the rate=0.2 and the third row count=5. On the fourth
##' ## row starts the next sequence with the values in the second row
##' ## in the edges data.frame. The stop value in the first column is
##' ## on row 10. As can be seen in column 6, there are no edge
##' ## properties for node=6.
##' ##        [,1] [,2]  [,3] [,4] [,5] [,6]
##' ##  [1,]  1.00    0  3.00  0.0  0.0   -1
##' ##  [2,]  0.20    1  0.20  0.2  0.2  NaN
##' ##  [3,]  5.00   50 10.00  5.0  5.0  NaN
##' ##  [4,]  2.00   -1  4.00  2.0  2.0  NaN
##' ##  [5,]  0.01  NaN  0.05  0.8  0.8  NaN
##' ##  [6,]  5.00  NaN 10.00  5.0  5.0  NaN
##' ##  [7,]  3.00  NaN -1.00 -1.0 -1.0  NaN
##' ##  [8,]  0.79  NaN   NaN  NaN  NaN  NaN
##' ##  [9,]  5.00  NaN   NaN  NaN  NaN  NaN
##' ## [10,] -1.00  NaN   NaN  NaN  NaN  NaN
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
    for (k in seq_len(length(j))) {
        m[cbind(i + k, edges$to)] <- edges[, j[k], drop = TRUE]
    }

    m
}
