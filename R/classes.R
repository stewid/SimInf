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

##' Class \code{SimInf_events}
##'
##' Class to hold data for scheduled events that modify the discrete
##' state of individuals in a node at a pre-defined time \code{t}.
##'
##' @slot E The \strong{select matrix} (sparse matrix of class
##'     \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}).
##'     Each row corresponds to a model compartment.
##'     \itemize{
##'       \item \strong{Sampling (Exit, Internal/External Transfer):}
##'         Non-zero entries in a column indicate which compartments
##'         individuals are sampled from. The values in \code{E[, select]}
##'         act as \strong{weights} for sampling individuals without
##'         replacement (probability proportional to weight).
##'       \item \strong{Targeting (Enter):}
##'         Non-zero entries in a column indicate which compartments
##'         new individuals are added to. The values in \code{E[, select]}
##'         act as \strong{weights} for distributing new individuals
##'         among the target compartments.
##'     }
##' @slot N The \strong{shift matrix} (integer matrix).  Determines
##'     how individuals are moved between compartments during
##'     \emph{enter}, \emph{internal transfer}, and \emph{external
##'     transfer} events.
##'     \itemize{
##'       \item Each row corresponds to a source compartment.
##'       \item Each column corresponds to a specific \code{shift}
##'       value.
##'       \item If \code{q <- shift}, the entry \code{N[p, q]} defines
##'         the \strong{offset} (number of rows to move) for
##'         individuals sampled from compartment \code{p}.
##'       \item The destination compartment is calculated as:
##'         \code{destination = p + N[p, q]}.
##'       \item Constraint: \code{1 <= destination <= number of
##'       compartments}.
##'     }
##' @slot event Integer vector specifying the event type for each row:
##'     \itemize{
##'       \item \code{0}: \emph{exit} (remove individuals).
##'       \item \code{1}: \emph{enter} (add individuals).
##'       \item \code{2}: \emph{internal transfer} (move within node).
##'       \item \code{3}: \emph{external transfer} (move between
##'       nodes).
##'     }
##'     Other values are reserved for future use.
##' @slot time Integer vector specifying the time step when each event
##'     occurs.
##' @slot node Integer vector specifying the \strong{source node} for
##'     the event.  For \emph{external transfer}, this is the node
##'     individuals are moved \emph{from}.  Range: \code{1 <= node <=
##'     number of nodes}.
##' @slot dest Integer vector specifying the \strong{destination node}
##'     for \emph{external transfer} events (individuals moved
##'     \emph{to}).  For other event types, this value is ignored
##'     (typically set to 0).  Range: \code{1 <= dest <= number of
##'     nodes}.
##' @slot n Integer vector specifying the \strong{number of
##'     individuals} affected by the event. Must be \code{n >= 0}.
##' @slot proportion Numeric vector. If \code{n[i] == 0}, the number
##'     of individuals is sampled from a binomial distribution using
##'     \code{proportion[i]} and the current population size in the
##'     selected compartments. Range: \code{0 <= proportion <= 1}.
##' @slot select Integer vector specifying which \strong{column} of
##'     the matrix \code{E} to use for sampling/targeting for each
##'     event.  The specific individuals are chosen based on the
##'     non-zero entries in \code{E[, select[i]]}.
##' @slot shift Integer vector specifying which \strong{column} of the
##'     matrix \code{N} to use for shifting individuals for each
##'     event.  Unused for \emph{exit} events.
##' @export
setClass(
    "SimInf_events",
    slots = c(E          = "dgCMatrix",
              N          = "matrix",
              event      = "integer",
              time       = "integer",
              node       = "integer",
              dest       = "integer",
              n          = "integer",
              proportion = "numeric",
              select     = "integer",
              shift      = "integer")
)

##' Class \code{"SimInf_model"}
##'
##' Class to handle data for the \code{SimInf_model}.
##' @template G-slot
##' @template S-slot
##' @template U-slot
##' @template U_sparse-slot
##' @slot V The result matrix for the real-valued continuous
##'     state. \code{V[, j]} contains the real-valued state of the
##'     system at \code{tspan[j]}. Numeric matrix
##'     (\eqn{N_n}\code{dim(ldata)[1]} \eqn{\times}
##'     \code{length(tspan)}).
##' @slot V_sparse If the model was configured to write the solution
##'     to a sparse matrix
##'     (\code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}) the
##'     \code{V_sparse} contains the data and \code{V} is empty. The
##'     layout of the data in \code{V_sparse} is identical to
##'     \code{V}.
##' @template ldata-slot
##' @template gdata-slot
##' @template tspan-slot
##' @template u0-slot
##' @slot v0 The initial value for the real-valued continuous state.
##'     Numeric matrix (\code{dim(ldata)[1]} \eqn{\times N_n}).
##' @slot events Scheduled events \code{\linkS4class{SimInf_events}}
##' @slot replicates Number of replicates of the model.
##' @slot C_code Character vector with optional model C code. If
##'     non-empty, the C code is written to a temporary C-file when
##'     the \code{run} method is called.  The temporary C-file is
##'     compiled and the resulting DLL is dynamically loaded.
##' @export
setClass(
    "SimInf_model",
    slots = c(G          = "dgCMatrix",
              S          = "dgCMatrix",
              U          = "matrix",
              U_sparse   = "dgCMatrix",
              ldata      = "matrix",
              gdata      = "numeric",
              tspan      = "numeric",
              u0         = "matrix",
              V          = "matrix",
              V_sparse   = "dgCMatrix",
              v0         = "matrix",
              events     = "SimInf_events",
              replicates = "integer",
              C_code     = "character")
)
