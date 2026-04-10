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
##'
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
##'
##' @slot event Integer vector specifying the event type for each row:
##'     \itemize{
##'       \item \code{0}: \emph{exit} (remove individuals).
##'       \item \code{1}: \emph{enter} (add individuals).
##'       \item \code{2}: \emph{internal transfer} (move within node).
##'       \item \code{3}: \emph{external transfer} (move between
##'       nodes).
##'     }
##'     Other values are reserved for future use.
##'
##' @slot time Integer vector specifying the time step when each event
##'     occurs.
##'
##' @slot node Integer vector specifying the \strong{source node} for
##'     the event.  For \emph{external transfer}, this is the node
##'     individuals are moved \emph{from}.  Range: \code{1 <= node <=
##'     number of nodes}.
##'
##' @slot dest Integer vector specifying the \strong{destination node}
##'     for \emph{external transfer} events (individuals moved
##'     \emph{to}).  For other event types, this value is ignored
##'     (typically set to 0).  Range: \code{1 <= dest <= number of
##'     nodes}.
##'
##' @slot n Integer vector specifying the \strong{number of
##'     individuals} affected by the event. Must be \code{n >= 0}.
##'
##' @slot proportion Numeric vector. If \code{n[i] == 0}, the number
##'     of individuals is sampled from a binomial distribution using
##'     \code{proportion[i]} and the current population size in the
##'     selected compartments. Range: \code{0 <= proportion <= 1}.
##'
##' @slot select Integer vector specifying which \strong{column} of
##'     the matrix \code{E} to use for sampling/targeting for each
##'     event.  The specific individuals are chosen based on the
##'     non-zero entries in \code{E[, select[i]]}.
##'
##' @slot shift Integer vector specifying which \strong{column} of the
##'     matrix \code{N} to use for shifting individuals for each
##'     event.  Unused for \emph{exit} events.
##'
##' @seealso \code{\link{SimInf_model}} for the main model class that
##'     holds the events.  \code{\link{SIR}}, \code{\link{SEIR}},
##'     \code{\link{SIS}}, \code{\link{SISe}} for examples of how
##'     events are passed to model constructors.
##'     \code{\link{SimInf_model}} (constructor) for details on the
##'     \code{E} and \code{N} matrices used to define event behavior.
##'     \code{\link{run}} for executing the simulation with scheduled
##'     events.  Vignette \code{"Scheduled events"} for detailed
##'     examples of defining and using events.
##'
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

##' Class \code{SimInf_model}
##'
##' The core class for storing the state, parameters, and results of a
##' stochastic simulation in \pkg{SimInf}. This class holds the model
##' definition (transition graphs, matrices), initial conditions,
##' scheduled events, and the simulation output.
##'
##' @slot G \strong{Dependency Graph} (sparse matrix, class
##'     \code{dgCMatrix}).  Indicates which transition rates need to
##'     be updated after a state transition occurs. A non-zero entry
##'     \code{G[i, j]} means that transition rate \code{i} must be
##'     recalculated if transition \code{j} occurs. This optimizes
##'     performance by avoiding unnecessary updates.  Dimensions:
##'     \eqn{N_t \times N_t}, where \eqn{N_t} is the number of
##'     transitions.
##'
##' @slot S \strong{State Transition Matrix} (sparse matrix, class
##'     \code{dgCMatrix}).  Defines the change in the state vector for
##'     each transition.  Executing transition \code{j} adds the
##'     column \code{S[, j]} to the state vector of the affected node.
##'     Dimensions: \eqn{N_c \times N_t}, where \eqn{N_c} is the
##'     number of compartments.
##'
##' @slot U \strong{Discrete State Result Matrix} (integer matrix).
##'     Contains the number of individuals in each compartment for
##'     every node at each time point in \code{tspan}.
##'     \itemize{
##'       \item \code{U[, j]}: State at time \code{tspan[j]}.
##'       \item Rows are ordered by node, then by compartment:
##'         \itemize{
##'           \item Rows \code{1:Nc}: Node 1.
##'           \item Rows \code{(Nc+1):(2*Nc)}: Node 2.
##'           \item ... and so on.
##'         }
##'     }
##'     Dimensions: \eqn{N_n \times N_c \times \text{length(tspan)}}.
##'     \emph{Note: If the model was run with sparse output, this slot
##'     is empty.}
##'
##' @slot U_sparse \strong{Sparse Discrete State Result} (sparse
##'     matrix, class \code{dgCMatrix}).  Contains the simulation
##'     results if the model was configured for sparse output.  The
##'     layout is identical to \code{U}, but stored as a sparse matrix
##'     to save memory.  \emph{Note: Only one of \code{U} or
##'     \code{U_sparse} will contain data.}
##'
##' @slot V \strong{Continuous State Result Matrix} (numeric matrix).
##'     Contains the values of continuous state variables (e.g.,
##'     environmental pathogen load) for every node at each time
##'     point.  Dimensions: \eqn{N_n \times N_{ld} \times
##'     \text{length(tspan)}}, where \eqn{N_{ld}} is the number of
##'     local data variables (continuous states).  \emph{Note: If
##'     sparse output was used, this slot is empty.}
##'
##' @slot V_sparse \strong{Sparse Continuous State Result} (sparse
##'     matrix, class \code{dgCMatrix}).  Contains the continuous
##'     state results if sparse output was enabled.  Layout identical
##'     to \code{V}.  \emph{Note: Only one of \code{V} or
##'     \code{V_sparse} will contain data.}
##'
##' @slot ldata \strong{Local Data Matrix} (numeric matrix).
##'     Parameters specific to each node (e.g., node-specific
##'     transmission rates).  Column \code{ldata[, j]} contains the
##'     local data vector for node \code{j}.  Passed to transition
##'     rate functions and post-step functions.  Dimensions:
##'     \eqn{N_{ld} \times N_n}.
##'
##' @slot gdata \strong{Global Data Vector} (numeric vector).
##'     Parameters common to all nodes (e.g., global recovery rate).
##'     Passed to transition rate functions and post-step functions.
##'
##' @slot tspan \strong{Time Span} (numeric vector).  Increasing time
##'     points where the state of each node is recorded.
##'
##' @slot u0 \strong{Initial State Matrix} (integer matrix).  Initial
##'     number of individuals in each compartment for every node.
##'     Dimensions: \eqn{N_c \times N_n}.
##'
##' @slot v0 \strong{Initial Continuous State Matrix} (numeric
##'     matrix).  Initial values for continuous state variables for
##'     every node.  Dimensions: \eqn{N_{ld} \times N_n}.
##'
##' @slot events \strong{Scheduled Events}
##'     (\code{\linkS4class{SimInf_events}}).  Object containing the
##'     schedule of discrete events (e.g., movements, births).
##'
##' @slot replicates \strong{Number of Replicates} (integer).  Number
##'     of parallel replicates simulated for this model (used in
##'     filtering algorithms).
##'
##' @slot C_code \strong{C Source Code} (character vector).  Optional
##'     C code defining the model's transition rates. If non-empty,
##'     this code is written to a temporary file, compiled, and loaded
##'     when \code{run()} is called.  Typically generated by
##'     \code{\link{mparse}}.
##'
##' @seealso \code{\link{SimInf_model}} (constructor) for creating
##'     model objects.  \code{\link{SIR}}, \code{\link{SEIR}},
##'     \code{\link{SIS}}, \code{\link{SISe}} for compartment model
##'     constructors that automatically set up the required slots.
##'     \code{\link{mparse}} for creating custom models using a simple
##'     string syntax.  \code{\link{run}} for executing the
##'     simulation.  \code{\link{trajectory}},
##'     \code{\link{prevalence}}, and \code{\link{plot}} for
##'     extracting and visualizing results.
##'     \code{\linkS4class{SimInf_events}} for details on the event
##'     schedule structure.
##'
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
