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

##' Class \code{"SimInf_events"}
##'
##' Class to hold data for scheduled events to modify the discrete
##' state of individuals in a node at a pre-defined time t.
##' @slot E Each row corresponds to one compartment in the model. The
##'     non-zero entries in a column indicate the compartments to
##'     include in an event. For the \emph{exit}, \emph{internal
##'     transfer} and \emph{external transfer} events, a non-zero
##'     entry indicates the compartments to sample individuals from,
##'     where the values in \code{E[, select]} are used as weights.
##'     Individuals are sampled without replacement with probability
##'     proportional to the weight in \code{E[, select]}. For the
##'     \emph{enter} event, the values in \code{E[, select]} are used
##'     as weights when determining which compartment to add
##'     individuals to. If the column \code{E[, select]} contains
##'     several non-zero entries, the compartment is sampled with
##'     probability proportional to the weight in \code{E[, select]}.
##'     \code{E} is sparse matrix of class
##'     \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}.
##' @slot N Determines how individuals in \emph{enter}, \emph{internal
##'     transfer} and \emph{external transfer} events are shifted to
##'     enter another compartment. Each row corresponds to one
##'     compartment in the model. The values in a column define how to
##'     move sampled individuals before adding them to the
##'     destination.  Let \code{q <- shift}, then each non-zero entry
##'     in \code{N[, q]} defines the number of rows to move sampled
##'     individuals from that compartment i.e., sampled individuals
##'     from compartment \code{p} are moved to compartment \code{N[p,
##'     q] + p}, where \code{1 <= N[p, q] + p <=
##'     N_compartments}. Which column to use for each event is
##'     specified by the \code{shift} vector (see below). \code{N} is
##'     an integer matrix.
##' @slot event Type of event: 0) \emph{exit}, 1) \emph{enter}, 2)
##'     \emph{internal transfer}, and 3) \emph{external transfer}.
##'     Other values are reserved for future event types and not
##'     supported by the current solvers. Integer vector.
##' @slot time Time of when the event occurs i.e., the event is
##'     processed when time is reached in the simulation.  \code{time}
##'     is an integer vector.
##' @slot node The node that the event operates on. Also the source
##'     node for an \emph{external transfer} event.  Integer vector.
##'     1 <= \code{node[i]} <= Number of nodes.
##' @slot dest The destination node for an \emph{external transfer}
##'     event i.e., individuals are moved from \code{node} to
##'     \code{dest}, where 1 <= \code{dest[i]} <= Number of nodes.
##'     Set \code{event = 0} for the other event types.  \code{dest}
##'     is an integer vector.
##' @slot n The number of individuals affected by the event. Integer
##'     vector.  n[i] >= 0.
##' @slot proportion If \code{n[i]} equals zero, the number of
##'     individuals affected by \code{event[i]} is calculated by
##'     sampling the number of individuals from a binomial
##'     distribution using the \code{proportion[i]} and the number of
##'     individuals in the compartments. Numeric vector.  0 <=
##'     proportion[i] <= 1.
##' @slot select To process \code{event[i]}, the compartments affected
##'     by the event are specified with \code{select[i]} together with
##'     the matrix \code{E}, where \code{select[i]} determines which
##'     column in \code{E} to use.  The specific individuals affected
##'     by the event are proportionally sampled from the compartments
##'     corresponding to the non-zero entries in the specified column
##'     in \code{E[, select[i]]}, where \code{select} is an integer
##'     vector.
##' @slot shift Determines how individuals in \emph{enter},
##'     \emph{internal transfer} and \emph{external transfer} events
##'     are shifted to enter another compartment. The sampled
##'     individuals are shifted according to column \code{shift[i]} in
##'     matrix \code{N} i.e., \code{N[, shift[i]]}, where \code{shift}
##'     is an integer vector.  Unused for \emph{exit} events.
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
