## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
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

##' Create a \code{SimInf_raster_model}
##'
##' Construct a low-level \code{SimInf_raster_model} object.  It is a
##' model where the nodes are not fixed at one position but can move
##' between cells on a raster.  This function is typically used
##' internally by model constructors (e.g., \code{mparse()}) or for
##' advanced usage where custom model definitions (e.g., user-provided
##' C code or non-standard matrices) are required.
##'
##' @param raster integer matrix with the landcover class in each
##'     cell.  Dimensions should be \code{nrow} \eqn{\times}
##'     \code{ncol} matching the spatial extent of the
##'     model. Landcover classes are zero-based integers (i.e., the
##'     first class is 0, second is 1, etc.). Row indices start at 1
##'     at the top and increase toward the bottom; column indices
##'     start at 1 at the left and increase to the right. Values must
##'     be non-negative integers; NA values are not supported.
##'
##' @param tr_type FIXME.
##'
##' @param G \strong{Dependency Graph}.  Indicates which transition
##'     rates need updating after a state transition.  Can be provided
##'     as a sparse matrix (class \code{dgCMatrix}) or a dense matrix.
##'     If a dense matrix is provided, it is automatically converted
##'     to a sparse format internally.  See
##'     \code{\linkS4class{SimInf_model}} for detailed matrix layout.
##'
##' @param S \strong{State Transition Matrix}.  Defines the change in
##'     the state vector for each transition.  Can be provided as a
##'     sparse matrix (class \code{dgCMatrix}) or a dense matrix.  If
##'     a dense matrix is provided, it is automatically converted to a
##'     sparse format internally.  See
##'     \code{\linkS4class{SimInf_model}} for detailed matrix layout.
##'
##' @param cell_S FIXME.
##'
##' @param tspan \strong{Time Span} (numeric or Date vector).
##'     Increasing time points for output. If \code{Date}, converted
##'     to days with names, where \code{tspan[1]} becomes the day of
##'     the year of the first year of \code{tspan}. The dates are
##'     added as names to the numeric vector.
##'
##' @param ldata \strong{Local Data}.
##'     Parameters specific to each node. Can be:
##'     \itemize{
##'       \item A \code{data.frame} with one row per node.
##'       \item A matrix where each column \code{ldata[, j]} is the
##'       data vector for node \code{j}.
##'     }
##'     Passed to transition rate and post-step functions.
##'
##' @param gdata \strong{Global Data} (numeric vector).  Parameters
##'     common to all nodes. Passed to transition rate and post-step
##'     functions.
##'
##' @param u0 \strong{Initial State}.  Initial number of individuals
##'     per compartment/node. Can be:
##'     \itemize{
##'       \item A matrix (\eqn{N_c \times N_n}).
##'       \item A \code{data.frame} with columns corresponding to
##'       compartments.
##'       \item Any object coercible to a \code{data.frame} (e.g., a
##'         named numeric vector will be coerced to a one-row
##'         \code{data.frame}).
##'     }
##'
##' @param v0 \strong{Initial Continuous State} (numeric matrix).
##'     Initial values for continuous states per node.
##'
##' @param C_code \strong{C Source Code} (character vector).  Optional
##'     C code for custom transition rates. If provided, it is
##'     compiled and loaded when \code{run()} is called.
##'
##' @export
SimInf_raster_model <- function(raster,
                                tr_type,
                                G,
                                S,
                                cell_S,
                                tspan,
                                ldata  = NULL,
                                gdata  = NULL,
                                u0     = NULL,
                                v0     = NULL,
                                C_code = NULL) {

    model <- SimInf_model(G      = G,
                          S      = S,
                          tspan  = tspan,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0,
                          C_code = C_code)

    model <- methods::as(model, "SimInf_raster_model")
    model@cell_S <- init_sparse_matrix(cell_S)
    model@raster <- init_x0(raster)
    model@tr_type <- as.integer(tr_type)

    model
}
