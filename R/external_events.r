## siminf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015  Stefan Engblom
## Copyright (C) 2015  Stefan Widgren
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

##' Check if wholenumbers
##'
##' Check that all values are wholenumbers, see example in integer {base}
##' @param x Value to check
##' @param tol Tolerance of the check
##' @return logical vector
##' @keywords internal
is_wholenumber <- function(x, tol = .Machine$double.eps^0.5)
{
    abs(x - round(x)) < tol
}

##' Class \code{"external_events"}
##'
##' Class to handle the external events
##' @section Slots:
##' \describe{
##'   \item{E}{
##'     Sparse matrix (\eqn{Ncompartments \times (4 * Nobs)}) of
##'     object class \code{"\linkS4class{dgCMatrix}"}. Each event
##'     type contains the number of observable states (Nobs) columns.
##'     Currently 4 types of events are implemented, see \code{event}.
##'     For the events EXIT_EVENT and EXTERNAL_TRANSFER_EVENT, a non-zero entry
##'     in element \code{i} of select column \code{j} indicate the
##'     compartments to sample individuals from.
##'     For the event ENTER_EVENT, all individuals enter first non-zero
##'     compartment, i.e. a non-zero entry in element \code{i} of select
##'     column \code{j}.
##'     For the INTERNAL_TRANSFER_EVENT, a non-zero entry
##'     in element \code{i} of select column \code{j} indicate the
##'     compartments to sample individuals from. The value of element
##'     \code{i} of select column \code{j} determines the target
##'     compartment. The target compartment for sampled individuals is
##'     given by adding the value.
##'   }
##'   \item{event}{
##'     Integer vector of length \code{len} with external events.
##'     The following four events are implemented; EXIT_EVENT = 0,
##'     ENTER_EVENT = 1, INTERNAL_TRANSFER_EVENT = 2, and
##'     EXTERNAL_TRANSFER_EVENT = 3.
##'   }
##'   \item{time}{
##'     Integer vector of length \code{len} with the time for external event.
##'   }
##'   \item{select}{
##'     Integer vector of length \code{len}. Column j in the event matrix
##'     that determines the hidden states to sample from.
##'   }
##'   \item{node}{
##'     Integer vector of length \code{len}. The source herd of the event i.
##'   }
##'   \item{dest}{
##'     Integer vector of length \code{len}. The dest herd of the event i.
##'   }
##'   \item{n}{
##'     Integer vector of length \code{len}. The number of individuals in the
##'     external event. n[i] >= 0.
##'   }
##'   \item{proportion}{
##'     Numeric vector of length \code{len}. If n[i] equals zero,
##'     then the number of individuals to sample is calculated by summing
##'     the number of individuals in the hidden states determined by
##'     select[i] and multiplying with the proportion.
##'     0 <= proportion[i] <= 1.
##'   }
##'   \item{len}{
##'     Number of scheduled external events.
##'   }
##' }
##' @name external_events-class
##' @docType class
##' @keywords classes
##' @keywords methods
##' @import methods
##' @import Matrix
##' @export
setClass("external_events",
         slots = c(E          = "dgCMatrix",
                   event      = "integer",
                   time       = "integer",
                   select     = "integer",
                   node       = "integer",
                   dest       = "integer",
                   n          = "integer",
                   proportion = "numeric",
                   len    = "integer"),
         prototype = list(len = 0L),
         validity = function(object) {
             errors <- character()

             if (!identical(length(unique(c(length(object@event),
                                            length(object@time),
                                            length(object@select),
                                            length(object@node),
                                            length(object@dest),
                                            length(object@n),
                                            length(object@proportion)))) , 1L)) {
                 errors <- c(errors, "All external events must have equal length.")
             }

             if (!identical(length(object@len), 1L)) {
                 errors <- c(errors, "Length of len must be equal to one.")
             }

             if (!identical(length(object@event),
                           object@len)) {
                 errors <- c(errors, "Length of external events must be equal to len.")
             }

             if (!all(object@time > 0)) {
                 errors <- c(errors, "time must be greater than 0")
             }

             if (any(object@event < 0, object@event > 3)) {
                 errors <- c(errors, "event must be in the range 0 <= event <= 3")
             }

             if (any(object@select < 0,
                     object@select >= (dim(object@E)[2] / 4))) {
                 errors <- c(errors, "select must be in the range 0 <= select < Nselect")
             }

             if (any(object@proportion < 0, object@proportion > 1)) {
                 errors <- c(errors, "prop must be in the range 0 <= prop <= 1")
             }

             if (length(errors) == 0) TRUE else errors
         }
)

##' Create S4 class \code{external_events}
##'
##' The argument events must be a \code{data.frame} with the following
##' columns:
##' \itemize{
##'   \item{event}{
##'     The event type
##'   }
##'   \item{time}{
##'     The time of the event
##'   }
##'   \item{select}{
##'     The internal category of the event
##'   }
##'   \item{node}{
##'     The herd of the event
##'   }
##'   \item{dest}{
##'     The destination herd of an animal movement
##'   }
##'   \item{n}{
##'     The number of individuals
##'   }
##'   \item{prop}{
##'     The proportion of individuals. If n[i] equals zero,
##'     then the number of individuals to sample is calculated by summing
##'     the number of individuals in the hidden states determined by
##'     select[i] and multiplying with the proportion. 0 <= p[i] <= 1.
##'   }
##' }
##' @param E Sparse matrix (\eqn{Ncompartments \times (4 * Nselect)}) of
##'        object class \code{"\linkS4class{dgCMatrix}"}.
##' @param events A \code{data.frame} with events.
##' @return S4 class \code{external_events}
##' @export
external_events <- function(E      = NULL,
                            events = NULL)
{
    ## Check E
    if (is.null(E)) {
        if (!is.null(events))
            stop("events is not NULL when E is NULL")
        E <- new("dgCMatrix")
    }

    ## Check events
    if (is.null(events)) {
        events <- data.frame(event      = as.integer(),
                             time       = as.integer(),
                             select     = as.integer(),
                             node       = as.integer(),
                             dest       = as.integer(),
                             n          = as.integer(),
                             proportion = as.numeric())
    }
    if (!is.data.frame(events))
        stop("events must be a data.frame")
    if (!identical(ncol(events), 7L))
        stop("Wrong dimension of events")
    if (!all(c("event", "time", "select", "node", "dest", "n", "proportion") %in% names(events)))
        stop("Missing columns in events")
    if (!is.numeric(events$event))
        stop("Columns in events must be numeric")
    if (!is.numeric(events$time))
        stop("Columns in events must be numeric")
    if (!is.numeric(events$select))
        stop("Columns in events must be numeric")
    if (!is.numeric(events$node))
        stop("Columns in events must be numeric")
    if (!is.numeric(events$dest))
        stop("Columns in events must be numeric")
    if (!is.numeric(events$n))
        stop("Columns in events must be numeric")
    if (!is.numeric(events$proportion))
        stop("Columns in events must be numeric")

    if (nrow(events)) {
        if (!all(is_wholenumber(events$event)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$time)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$select)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$node)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$dest)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$n)))
            stop("Columns in events must be integer")
    }

    events <- events[order(events$time, events$event, events$node),]

    return(new("external_events",
               E          = E,
               event      = as.integer(events$event),
               time       = as.integer(events$time),
               select     = as.integer(events$select),
               node       = as.integer(events$node),
               dest       = as.integer(events$dest),
               n          = as.integer(events$n),
               proportion = as.numeric(events$proportion),
               len        = nrow(events)))
}

##' Brief summary of \code{external_events}
##'
##' @aliases show,external_events-methods
##' @docType methods
##' @param object The external_events \code{object}
##' @return None (invisible 'NULL').
##' @keywords methods
##' @export
setMethod("show",
          signature(object = "external_events"),
          function (object)
          {
              cat("\nExternal events:\n")
              cat(sprintf("E: %i x %i\n", dim(object@E)[1], dim(object@E)[2]))

              if (length(object@event)) {
                  cat(sprintf("event: 1 x %i\n", length(object@event)))
              } else {
                  cat("event: 0 x 0\n")
              }

              if (length(object@time)) {
                  cat(sprintf("time: 1 x %i\n", length(object@time)))
              } else {
                  cat("time: 0 x 0\n")
              }

              if (length(object@select)) {
                  cat(sprintf("select: 1 x %i\n", length(object@select)))
              } else {
                  cat("select: 0 x 0\n")
              }

              if (length(object@node)) {
                  cat(sprintf("node: 1 x %i\n", length(object@node)))
              } else {
                  cat("node: 0 x 0\n")
              }

              if (length(object@dest)) {
                  cat(sprintf("dest: 1 x %i\n", length(object@dest)))
              } else {
                  cat("dest: 0 x 0\n")
              }

              if (length(object@n)) {
                  cat(sprintf("n: 1 x %i\n", length(object@n)))
              } else {
                  cat("n: 0 x 0\n")
              }

              if (length(object@proportion)) {
                  cat(sprintf("proportion: 1 x %i\n", length(object@proportion)))
              } else {
                  cat("proportion: 0 x 0\n")
              }

              if (length(object@len)) {
                  cat(sprintf("len: %i\n", object@len[1]))
              } else {
                  cat("len: 0\n")
              }
          }
)
