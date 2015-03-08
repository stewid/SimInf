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
##'     Currently 4 types of events are implemented, see \code{ext_event}.
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
##'   \item{ext_event}{
##'     Integer vector of length ext_len with external events.
##'     The following four events are implemented; EXIT_EVENT = 0,
##'     ENTER_EVENT = 1, INTERNAL_TRANSFER_EVENT = 2, and
##'     EXTERNAL_TRANSFER_EVENT = 3.
##'   }
##'   \item{ext_time}{
##'     Integer vector of length ext_len with the time for external event.
##'   }
##'   \item{ext_select}{
##'     Integer vector of length ext_len. Column j in the event matrix
##'     that determines the hidden states to sample from.
##'   }
##'   \item{ext_node}{
##'     Integer vector of length ext_len. The source herd of the event i.
##'   }
##'   \item{ext_dest}{
##'     Integer vector of length ext_len. The dest herd of the event i.
##'   }
##'   \item{ext_n}{
##'     Integer vector of length ext_len. The number of individuals in the external event. ext_n[i] >= 0.
##'   }
##'   \item{ext_p}{
##'     Numeric vector of length ext_len. If ext_n[i] equals zero,
##'     then the number of individuals to sample is calculated by summing
##'     the number of individuals in the hidden states determined by
##'     ext_select[i] and multiplying with the proportion. 0 <= ext_p[i] <= 1.
##'   }
##'   \item{ext_len}{
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
                   ext_event  = "integer",
                   ext_time   = "integer",
                   ext_select = "integer",
                   ext_node   = "integer",
                   ext_dest   = "integer",
                   ext_n      = "integer",
                   ext_p      = "numeric",
                   ext_len    = "integer"),
         prototype = list(ext_len = 0L),
         validity = function(object) {
             errors <- character()

             if (!identical(length(unique(c(length(object@ext_event),
                                            length(object@ext_time),
                                            length(object@ext_select),
                                            length(object@ext_node),
                                            length(object@ext_dest),
                                            length(object@ext_n),
                                            length(object@ext_p)))) , 1L)) {
                 errors <- c(errors, "All external events must have equal length.")
             }

             if (!identical(length(object@ext_len), 1L)) {
                 errors <- c(errors, "Length of ext_len must be equal to one.")
             }

             if (!identical(length(object@ext_event),
                           object@ext_len)) {
                 errors <- c(errors, "Length of external events must be equal to ext_len.")
             }

             if (!all(object@ext_time > 0)) {
                 errors <- c(errors, "time must be greater than 0")
             }

             if (any(object@ext_event < 0, object@ext_event > 3)) {
                 errors <- c(errors, "event must be in the range 0 <= event <= 3")
             }

             if (any(object@ext_select < 0,
                     object@ext_select >= (dim(object@E)[2] / 4))) {
                 errors <- c(errors, "select must be in the range 0 <= select < Nselect")
             }

             if (any(object@ext_p < 0, object@ext_p > 1)) {
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
external_events <- function(E       = NULL,
                            events  = NULL)
{
    ## Check E
    if (is.null(E)) {
        if (!is.null(events))
            stop("events is not NULL when E is NULL")
        E <- new("dgCMatrix")
    }

    ## Check events
    if (is.null(events)) {
        events <- data.frame(event  = as.integer(),
                             time   = as.integer(),
                             select = as.integer(),
                             node   = as.integer(),
                             dest   = as.integer(),
                             n      = as.integer(),
                             prop   = as.numeric())
    }
    if (!is.data.frame(events))
        stop("events must be a data.frame")
    if (!identical(ncol(events), 7L))
        stop("Wrong dimension of events")
    if (!all(c("event", "time", "select", "node", "dest", "n", "prop") %in% names(events)))
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
    if (!is.numeric(events$prop))
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
               ext_event  = as.integer(events$event),
               ext_time   = as.integer(events$time),
               ext_select = as.integer(events$select),
               ext_node   = as.integer(events$node),
               ext_dest   = as.integer(events$dest),
               ext_n      = as.integer(events$n),
               ext_p      = as.numeric(events$prop),
               ext_len    = nrow(events)))
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

              if (length(object@ext_event)) {
                  cat(sprintf("ext_event: 1 x %i\n", length(object@ext_event)))
              } else {
                  cat("ext_event: 0 x 0\n")
              }

              if (length(object@ext_time)) {
                  cat(sprintf("ext_time: 1 x %i\n", length(object@ext_time)))
              } else {
                  cat("ext_time: 0 x 0\n")
              }

              if (length(object@ext_select)) {
                  cat(sprintf("ext_select: 1 x %i\n", length(object@ext_select)))
              } else {
                  cat("ext_select: 0 x 0\n")
              }

              if (length(object@ext_node)) {
                  cat(sprintf("ext_node: 1 x %i\n", length(object@ext_node)))
              } else {
                  cat("ext_node: 0 x 0\n")
              }

              if (length(object@ext_dest)) {
                  cat(sprintf("ext_dest: 1 x %i\n", length(object@ext_dest)))
              } else {
                  cat("ext_dest: 0 x 0\n")
              }

              if (length(object@ext_n)) {
                  cat(sprintf("ext_n: 1 x %i\n", length(object@ext_n)))
              } else {
                  cat("ext_n: 0 x 0\n")
              }

              if (length(object@ext_p)) {
                  cat(sprintf("ext_p: 1 x %i\n", length(object@ext_p)))
              } else {
                  cat("ext_p: 0 x 0\n")
              }

              if (length(object@ext_len)) {
                  cat(sprintf("ext_len: %i\n", object@ext_len[1]))
              } else {
                  cat("ext_len: 0\n")
              }
          }
)
