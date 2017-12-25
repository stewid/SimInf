## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
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

##' Check if wholenumbers
##'
##' Check that all values are wholenumbers, see example in integer {base}
##' @param x Value to check
##' @param tol Tolerance of the check
##' @return logical vector
##' @noRd
is_wholenumber <- function(x, tol = .Machine$double.eps^0.5)
{
    abs(x - round(x)) < tol
}

##' Class \code{"SimInf_events"}
##'
##' Class to handle the scheduled events
##' @section Slots:
##' \describe{
##'   \item{E}{
##'     Each row corresponds to one compartment in the model. The
##'     non-zero entries in a column indicates the compartments to
##'     include in an event.  For the \emph{exit}, \emph{internal
##'     transfer} and \emph{external transfer} events, a non-zero
##'     entry indicate the compartments to sample individuals from.
##'     For the \emph{enter} event, all individuals enter first
##'     non-zero compartment. Sparse matrix of object class
##'     \code{\linkS4class{dgCMatrix}}.
##'   }
##'   \item{N}{
##'     Each row represents one compartment in the model and the
##'     values determine how to move sampled individuals in
##'     \emph{internal transfer} and \emph{external transfer} events.
##'   }
##'   \item{event}{
##'     Type of event: 0) \emph{exit}, 1) \emph{enter}, 2)
##'     \emph{internal transfer}, and 3) \emph{external transfer}.
##'     Other values are reserved for future event types and not
##'     supported by the current default core solver. Integer vector.
##'   }
##'   \item{time}{
##'     Time for the event. Integer vector.
##'   }
##'   \item{node}{
##'     The node that the event operates on. Also the source node for
##'     an \emph{external transfer} event.  Integer vector.
##'     1 <= \code{node[i]} <= Number of nodes.
##'   }
##'   \item{dest}{
##'     The destination node for an \emph{external transfer} event.
##'     Should be \code{0} for the other event types.
##'     Integer vector. dest[i] >= 0.
##'   }
##'   \item{n}{
##'     The number of individuals affected by the event. Integer vector.
##'     n[i] >= 0.
##'   }
##'   \item{proportion}{
##'     If \code{n[i]} equals zero, the number of individuals affected by
##'     \code{event[i]} is calculated by summing the number of individuls
##      in the compartments determined by \code{select[i]} and multiplying
##'     with \code{proportion[i]}. Numeric vector.
##'     0 <= proportion[i] <= 1.
##'   }
##'   \item{select}{
##'     The column \code{j} in the event matrix \code{E} that
##'     determines the compartments that the event operates
##'     on. Integer vector.
##'   }
##'   \item{shift}{
##'     The column \code{k} in the shift matrix \code{N} that
##'     determines how individuals in \emph{internal transfer} and
##'     \emph{external transfer} events are shifted to enter another
##'     compartment.  Should be \code{0} for the other event types.
##'     Integer vector.
##'   }
##' }
##' @export
setClass("SimInf_events",
         slots = c(E          = "dgCMatrix",
                   N          = "matrix",
                   event      = "integer",
                   time       = "integer",
                   node       = "integer",
                   dest       = "integer",
                   n          = "integer",
                   proportion = "numeric",
                   select     = "integer",
                   shift      = "integer"),
         validity = function(object) {
             errors <- character()

             if (!identical(length(unique(c(length(object@event),
                                            length(object@time),
                                            length(object@node),
                                            length(object@dest),
                                            length(object@n),
                                            length(object@proportion),
                                            length(object@select),
                                            length(object@shift)))) , 1L)) {
                 errors <- c(errors,
                             "All scheduled events must have equal length.")
             }

             if (!all(object@time > 0)) {
                 errors <- c(errors,
                             "time must be greater than 0")
             }

             if (any(object@event < 0, object@event > 3)) {
                 errors <- c(errors,
                             "event must be in the range 0 <= event <= 3")
             }

             if (any(object@node < 1)) {
                 errors <- c(errors,
                             "'node' must be greater or equal to 1")
             }

             if (any(object@dest[object@event == 3] < 1)) {
                 errors <- c(errors,
                             "'dest' must be greater or equal to 1")
             }

             if (any(object@proportion < 0, object@proportion > 1)) {
                 errors <- c(errors,
                             "prop must be in the range 0 <= prop <= 1")
             }

             if (any(object@select < 1, object@select > dim(object@E)[2])) {
                 errors <- c(errors,
                             "select must be in the range 1 <= select <= Nselect")
             }

             if (any(object@shift[object@event == 2] < 1)) {
                 errors <- c(errors,
                             "'shift' must be greater or equal to 1")
             }

             if (length(errors) == 0) TRUE else errors
         }
)

##' Create a \code{\linkS4class{SimInf_events}} object
##'
##' The argument events must be a \code{data.frame} with the following
##' columns:
##' \describe{
##'   \item{event}{
##'     Type of event: 0) \emph{exit}, 1) \emph{enter}, 2)
##'     \emph{internal transfer}, and 3) \emph{external transfer}.
##'     Other values are reserved for future event types and not
##'     supported by the current default core solver.
##'   }
##'   \item{time}{
##'     Time for the event. Can be either an \code{integer} or a
##'     \code{Date} vector.  A \code{Date} vector is coerced to a
##'     numeric vector as days, where \code{t0} determines the offset
##'     to match the time of the events to the model \code{tspan}
##'     vector.
##'   }
##'   \item{node}{
##'     The node that the event operates on. Also the source node for
##'     an \emph{external transfer} event.
##'     1 <= \code{node[i]} <= Number of nodes.
##'   }
##'   \item{dest}{
##'     The destination node for an \emph{external transfer} event.
##'     Should be \code{0} for the other event types. dest[i] >= 0.
##'   }
##'   \item{n}{
##'     The number of individuals affected by the event. n[i] >= 0.
##'   }
##'   \item{proportion}{
##'     If \code{n[i]} equals zero, the number of individuals affected by
##'     \code{event[i]} is calculated by summing the number of individuls
##      in the compartments determined by \code{select[i]} and multiplying
##'     with \code{proportion[i]}. 0 <= proportion[i] <= 1.
##'   }
##'   \item{select}{
##'     The column \code{j} in the event matrix \code{E} that
##'     determines the compartments that the event operates on.
##'   }
##'   \item{shift}{
##'     The column \code{k} in the shift matrix \code{N} that
##'     determines how individuals in \emph{internal transfer} and
##'     \emph{external transfer} events are shifted to enter another
##'     compartment.  Should be \code{0} for the other event types.
##'   }
##' }
##'
##' @param E Sparse matrix of object class
##'     \code{\linkS4class{dgCMatrix}} that selects the states to
##'     include for sampling in an event.
##' @param N Sparse matrix of object class
##'     \code{\linkS4class{dgCMatrix}} that determines how to shift
##'     the states in an \code{INTERNAL_TRANSFER_EVENT}.
##' @param events A \code{data.frame} with events.
##' @param t0 If \code{events$time} is a \code{Date} vector, then
##'     \code{t0} determines the offset to match the time of the
##'     events to the model \code{tspan} vector, see details. If
##'     \code{events$time} is a numeric vector, then \code{t0} must be
##'     \code{NULL}.
##' @return S4 class \code{SimInf_events}
##' @export
##' @importFrom methods as
##' @importFrom methods is
##' @importFrom methods new
SimInf_events <- function(E      = NULL,
                          N      = NULL,
                          events = NULL,
                          t0     = NULL)
{
    ## Check E
    if (is.null(E)) {
        if (!is.null(events))
            stop("events is not NULL when E is NULL")
        E <- new("dgCMatrix")
    } else if (!is(E, "dgCMatrix")) {
        E <- as(E, "dgCMatrix")
    }

    ## Check N
    if (is.null(N)) {
        if (!is.null(events))
            stop("events is not NULL when N is NULL")
        N <- matrix(integer(0), nrow = 0, ncol = 0)
    }
    if (!all(is.matrix(N), is.numeric(N)))
        stop("'N' must be an integer matrix")
    if (!is.integer(N)) {
        if (!all(is_wholenumber(N)))
            stop("'N' must be an integer matrix")
        storage.mode(N) <- "integer"
    }

    ## Check events
    if (is.null(events)) {
        events <- data.frame(event      = as.integer(),
                             time       = as.integer(),
                             node       = as.integer(),
                             dest       = as.integer(),
                             n          = as.integer(),
                             proportion = as.numeric(),
                             select     = as.integer(),
                             shift      = as.integer())
    }
    if (!is.data.frame(events))
        stop("events must be a data.frame")
    if (!identical(ncol(events), 8L))
        stop("Wrong dimension of events")
    if (!all(c("event", "time", "node", "dest",
               "n", "proportion", "select",
               "shift") %in% names(events))) {
        stop("Missing columns in events")
    }

    ## Check time
    if (nrow(events)) {
        if (is(events$time, "Date")) {
            if (is.null(t0))
                stop("Missing 't0'")
            if (!all(identical(length(t0), 1L), is.numeric(t0)))
                stop("Invalid 't0'")
            events$time <- as.numeric(events$time) - t0
        } else if (!is.null(t0)) {
            stop("Invalid 't0'")
        }
    }

    if (!all(is.numeric(events$event), is.numeric(events$time),
             is.numeric(events$node), is.numeric(events$dest),
             is.numeric(events$n), is.numeric(events$proportion),
             is.numeric(events$select))) {
        stop("Columns in events must be numeric")
    }

    if (nrow(events)) {
        if (!all(is_wholenumber(events$event)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$time)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$node)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$dest)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$n)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$select)))
            stop("Columns in events must be integer")
        if (!all(is_wholenumber(events$shift)))
            stop("Columns in events must be integer")
    }

    events <- events[order(events$time, events$event, events$select),]

    new("SimInf_events",
        E          = E,
        N          = N,
        event      = as.integer(events$event),
        time       = as.integer(events$time),
        node       = as.integer(events$node),
        dest       = as.integer(events$dest),
        n          = as.integer(events$n),
        proportion = as.numeric(events$proportion),
        select     = as.integer(events$select),
        shift      = as.integer(events$shift))
}

setAs(from = "SimInf_events",
      to = "data.frame",
      def = function(from)
      {
          data.frame(event = from@event,
                     time = from@time,
                     node = from@node,
                     dest = from@dest,
                     n = from@n,
                     proportion = from@proportion,
                     select = from@select,
                     shift = from@shift)
      }
)

##' Plot scheduled events
##'
##' @param x the time points of the events.
##' @param y the number of events over time.
##' @param events the event type to plot.
##' @param frame.plot a logical indicating whether a box should be
##'     drawn around the plot.
##' @param ... additional arguments affecting the plot.
##' @importFrom graphics plot
##' @importFrom graphics mtext
##' @noRd
plot_SimInf_events <- function(x,
                                  y,
                                  events = c("Exit",
                                             "Enter",
                                             "Internal transfer",
                                             "External transfer"),
                                  frame.plot,
                                  ...)
{
    events <- match.arg(events)
    i <- switch(events,
                "Exit" = "0",
                "Enter" = "1",
                "Internal transfer" = "2",
                "External transfer" = "3")

    if (length(x)) {
        ylim <- c(0, max(y))

        if (i %in% rownames(y)) {
            y <- y[i, ]
        } else {
            y <- rep(0, length(x))
        }

        plot(x, y, type = "l", ylim = ylim, xlab = "",
             ylab = "", frame.plot = frame.plot, ...)
    } else {
        plot(0, 0, type = "n", xlab = "", ylab = "",
             frame.plot = frame.plot, ...)
    }

    mtext(events, side = 3, line = 0)
    mtext("Individuals", side = 2, line = 2)
    mtext("Time", side = 1, line = 2)
}

##' Display the distribution of scheduled events over time
##'
##' @param x The events data to plot.
##' @param frame.plot Draw a frame around each plot. Default is FALSE.
##' @param ... Additional arguments affecting the plot
##' @aliases plot,SimInf_events-method
##' @export
##' @importFrom graphics par
##' @importFrom stats xtabs
setMethod("plot",
          signature(x = "SimInf_events"),
          function(x, frame.plot = FALSE, ...)
          {
              savepar <- par(mfrow = c(2, 2),
                             oma = c(1, 1, 2, 0),
                             mar = c(4, 3, 1, 1))
              on.exit(par(savepar))

              yy <- xtabs(n ~ event + time,
                          cbind(event = x@event, time = x@time, n = x@n))
              xx <- as.integer(colnames(yy))

              ## Plot events
              plot_SimInf_events(xx, yy, "Exit", frame.plot, ...)
              plot_SimInf_events(xx, yy, "Enter", frame.plot, ...)
              plot_SimInf_events(xx, yy, "Internal transfer", frame.plot, ...)
              plot_SimInf_events(xx, yy, "External transfer", frame.plot, ...)
          }
)

##' Brief summary of \code{SimInf_events}
##'
##' Shows the number of scheduled events.
##' @param object The SimInf_events \code{object}
##' @return None (invisible 'NULL').
##' @export
##' @importFrom methods show
setMethod("show",
          signature(object = "SimInf_events"),
          function (object)
          {
              cat(sprintf("Number of scheduled events: %i\n",
                          length(object@event)))
          }
)

##' Detailed summary of a \code{SimInf_events} object
##'
##' Shows the number of scheduled events and the number of scheduled
##' events per event type.
##' @param object The \code{SimInf_events} object
##' @param ... Additional arguments affecting the summary produced.
##' @return None (invisible 'NULL').
##' @export
setMethod("summary",
          signature(object = "SimInf_events"),
          function(object, ...)
          {
              n <- length(object@event)
              cat(sprintf("Number of scheduled events: %i\n", n))

              ## Summarise exit events
              i <- which(object@event == 0L)
              if (length(i) > 0) {
                  cat(sprintf(" - Exit: %i (n: min = %i max = %i avg = %.1f)\n",
                              length(i),
                              min(object@n[i]),
                              max(object@n[i]),
                              mean(object@n[i])))
              } else {
                  cat(" - Exit: 0\n")
              }

              ## Summarise enter events
              i <- which(object@event == 1L)
              if (length(i) > 0) {
                  cat(sprintf(" - Enter: %i (n: min = %i max = %i avg = %.1f)\n",
                              length(i),
                              min(object@n[i]),
                              max(object@n[i]),
                              mean(object@n[i])))
              } else {
                  cat(" - Enter: 0\n")
              }

              ## Summarise internal transfer events
              i <- which(object@event == 2L)
              if (length(i) > 0) {
                  cat(sprintf(" - Internal transfer: %i (n: min = %i max = %i avg = %.1f)\n",
                              length(i),
                              min(object@n[i]),
                              max(object@n[i]),
                              mean(object@n[i])))
              } else {
                  cat(" - Internal transfer: 0\n")
              }

              ## Summarise external transfer events
              i <- which(object@event == 3L)
              if (length(i) > 0) {
                  cat(sprintf(" - External transfer: %i (n: min = %i max = %i avg = %.1f)\n",
                              length(i),
                              min(object@n[i]),
                              max(object@n[i]),
                              mean(object@n[i])))
              } else {
                  cat(" - External transfer: 0\n")
              }
          }
)
