## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2019 Stefan Widgren
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
##' @section Slots:
##' \describe{
##'   \item{E}{
##'     Each row corresponds to one compartment in the model. The
##'     non-zero entries in a column indicates the compartments to
##'     include in an event.  For the \emph{exit}, \emph{internal
##'     transfer} and \emph{external transfer} events, a non-zero
##'     entry indicate the compartments to sample individuals from.
##'     For the \emph{enter} event, all individuals enter first
##'     non-zero compartment. \code{E} is sparse matrix of class
##'     \code{\linkS4class{dgCMatrix}}.
##'   }
##'   \item{N}{
##'      Determines how individuals in \emph{internal transfer} and
##'      \emph{external transfer} events are shifted to enter another
##'      compartment.  Each row corresponds to one compartment in the
##'      model.  The values in a column are added to the current
##'      compartment of sampled individuals to specify the destination
##'      compartment, for example, a value of \code{1} in an entry
##'      means that sampled individuals in this compartment are moved
##'      to the next compartment.  Which column to use for each event
##'      is specified by the \code{shift} vector (see below).
##'      \code{N} is an integer matrix.
##'   }
##'   \item{event}{
##'     Type of event: 0) \emph{exit}, 1) \emph{enter}, 2)
##'     \emph{internal transfer}, and 3) \emph{external transfer}.
##'     Other values are reserved for future event types and not
##'     supported by the current solvers. Integer vector.
##'   }
##'   \item{time}{
##'     Time of when the event occurs i.e., the event is processed
##'     when time is reached in the simulation.  \code{time} is an
##'     integer vector.
##'   }
##'   \item{node}{
##'     The node that the event operates on. Also the source node for
##'     an \emph{external transfer} event.  Integer vector.
##'     1 <= \code{node[i]} <= Number of nodes.
##'   }
##'   \item{dest}{
##'     The destination node for an \emph{external transfer} event
##'     i.e., individuals are moved from \code{node} to \code{dest},
##'     where 1 <= \code{dest[i]} <= Number of nodes.  Set \code{event
##'     = 0} for the other event types.  \code{dest} is an integer
##'     vector.
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
##'     To process \code{event[i]}, the compartments affected by the
##'     event are specified with \code{select[i]} together with the
##'     matrix \code{E}, where \code{select[i]} determines which
##'     column in \code{E} to use.  The specific individuals affected
##'     by the event are proportionally sampled from the compartments
##'     corresponding to the non-zero entries in the specified column
##'     in \code{E[, select[i]]}, where \code{select} is an integer
##'     vector.
##'   }
##'   \item{shift}{
##'     Determines how individuals in \emph{internal transfer} and
##'     \emph{external transfer} events are shifted to enter another
##'     compartment.  The sampled individuals are shifted according to
##'     column \code{shift[i]} in matrix \code{N} i.e., \code{N[,
##'     shift[i]]}, where \code{shift} is an integer vector.  See
##'     above for a description of \code{N}. Unsued for the other
##'     event types.
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
                   shift      = "integer"))

## Check if the SimInf_events object is valid.
valid_SimInf_events_object <- function(object)
{
    ## Check that E and N have identical compartments
    if ((dim(object@E)[1] > 0) && (dim(object@N)[1] > 0)) {
        if (any(is.null(rownames(object@E)), is.null(rownames(object@N))))
            return("'E' and 'N' must have rownames matching the compartments.")
        if (!identical(rownames(object@E), rownames(object@N)))
            return("'E' and 'N' must have identical compartments.")
    }

    if (!identical(length(unique(c(length(object@event),
                                   length(object@time),
                                   length(object@node),
                                   length(object@dest),
                                   length(object@n),
                                   length(object@proportion),
                                   length(object@select),
                                   length(object@shift)))), 1L)) {
        return("All scheduled events must have equal length.")
    }

    if (!all(object@time > 0))
        return("time must be greater than 0.")

    if (any(object@event < 0, object@event > 3))
        return("event must be in the range 0 <= event <= 3.")

    if (any(object@node < 1))
        return("'node' must be greater or equal to 1.")

    if (any(object@dest[object@event == 3] < 1))
        return("'dest' must be greater or equal to 1.")

    if (any(object@proportion < 0, object@proportion > 1))
        return("prop must be in the range 0 <= prop <= 1.")

    if (any(object@select < 1, object@select > dim(object@E)[2]))
        return("select must be in the range 1 <= select <= Nselect.")

    if (any(object@shift[object@event == 2] < 1))
        return("'shift' must be greater or equal to 1.")

    TRUE
}

## Assign the function as the validity method for the class.
setValidity("SimInf_events", valid_SimInf_events_object)

##' Create a \code{\linkS4class{SimInf_events}} object
##'
##' The argument events must be a \code{data.frame} with the following
##' columns:
##' \describe{
##'   \item{event}{
##'     Four event types are supported by the current solvers:
##'     \emph{exit}, \emph{enter}, \emph{internal transfer}, and
##'     \emph{external transfer}.  When assigning the events, they can
##'     either be coded as a numerical value or a character string:
##'     \emph{exit;} \code{0} or \code{'exit'}, \emph{enter;} \code{1}
##'     or \code{'enter'}, \emph{internal transfer;} \code{2} or
##'     \code{'intTrans'}, and \emph{external transfer;} \code{3} or
##'     \code{'extTrans'}.  Internally in \pkg{SimInf}, the event type
##'     is coded as a numerical value.
##'   }
##'   \item{time}{
##'     When the event occurs i.e., the event is processed when time
##'     is reached in the simulation. Can be either an \code{integer}
##'     or a \code{Date} vector.  A \code{Date} vector is coerced to a
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
##'     The destination node for an \emph{external transfer} event
##'     i.e., individuals are moved from \code{node} to \code{dest},
##'     where 1 <= \code{dest[i]} <= Number of nodes.  Set \code{event
##'     = 0} for the other event types.  \code{dest} is an integer
##'     vector.
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
##'     To process \code{event[i]}, the compartments affected by the
##'     event are specified with \code{select[i]} together with the
##'     matrix \code{E}, where \code{select[i]} determines which
##'     column in \code{E} to use.  The specific individuals affected
##'     by the event are proportionally sampled from the compartments
##'     corresponding to the non-zero entries in the specified column
##'     in \code{E[, select[i]]}, where \code{select} is an integer
##'     vector.
##'   }
##'   \item{shift}{
##'     Determines how individuals in \emph{internal transfer} and
##'     \emph{external transfer} events are shifted to enter another
##'     compartment.  The sampled individuals are shifted according to
##'     column \code{shift[i]} in matrix \code{N} i.e., \code{N[,
##'     shift[i]]}, where \code{shift} is an integer vector.  See
##'     above for a description of \code{N}. Unsued for the other
##'     event types.
##'   }
##' }
##'
##' @param E Each row corresponds to one compartment in the model. The
##'     non-zero entries in a column indicates the compartments to
##'     include in an event.  For the \emph{exit}, \emph{internal
##'     transfer} and \emph{external transfer} events, a non-zero
##'     entry indicate the compartments to sample individuals from.
##'     For the \emph{enter} event, all individuals enter first
##'     non-zero compartment. \code{E} is sparse matrix of class
##'     \code{\linkS4class{dgCMatrix}}.
##' @param N Determines how individuals in \emph{internal transfer}
##'     and \emph{external transfer} events are shifted to enter
##'     another compartment.  Each row corresponds to one compartment
##'     in the model.  The values in a column are added to the current
##'     compartment of sampled individuals to specify the destination
##'     compartment, for example, a value of \code{1} in an entry
##'     means that sampled individuals in this compartment are moved
##'     to the next compartment.  Which column to use for each event
##'     is specified by the \code{shift} vector (see below).  \code{N}
##'     is an integer matrix.
##' @param events A \code{data.frame} with events.
##' @param t0 If \code{events$time} is a \code{Date} vector, then
##'     \code{t0} determines the offset to match the time of the
##'     events to the model \code{tspan} vector, see details. If
##'     \code{events$time} is a numeric vector, then \code{t0} must be
##'     \code{NULL}.
##' @return S4 class \code{SimInf_events}
##' @include check_arguments.R
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
            stop("events is not NULL when E is NULL.", call. = FALSE)
        E <- new("dgCMatrix")
    } else if (!is(E, "dgCMatrix")) {
        E <- as(E, "dgCMatrix")
    }

    ## Check N
    N <- check_N(N)

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
        stop("events must be a data.frame.", call. = FALSE)
    if (!all(c("event", "time", "node", "dest",
               "n", "proportion", "select",
               "shift") %in% names(events))) {
        stop("Missing columns in events.", call. = FALSE)
    }

    ## Do we have to recode the event type as a numerical value
    if (is.character(events$event) || is.factor(events$event)) {
        if (!all(events$event %in% c("enter", "exit", "extTrans", "intTrans"))) {
            stop("'event' type must be 'enter', 'exit', 'extTrans' or 'intTrans'.",
                 call. = FALSE)
        }

        ## Find indices to 'enter', 'internal transfer' and 'external
        ## transfer' events.
        i_enter <- which(events$event == "enter")
        i_intTrans <- which(events$event == "intTrans")
        i_extTrans <- which(events$event == "extTrans")

        ## Replace the character event type with a numerical value.
        events$event <- rep(0L, nrow(events))
        events$event[i_enter] <- 1L
        events$event[i_intTrans] <- 2L
        events$event[i_extTrans] <- 3L
    }

    ## Check time
    if (nrow(events)) {
        if (is(events$time, "Date")) {
            if (is.null(t0))
                stop("Missing 't0'.", call. = FALSE)
            if (!all(identical(length(t0), 1L), is.numeric(t0)))
                stop("Invalid 't0'.", call. = FALSE)
            events$time <- as.numeric(events$time) - t0
        } else if (!is.null(t0)) {
            stop("Invalid 't0'.", call. = FALSE)
        }
    }

    if (!all(is.numeric(events$event), is.numeric(events$time),
             is.numeric(events$node), is.numeric(events$dest),
             is.numeric(events$n), is.numeric(events$proportion),
             is.numeric(events$select))) {
        stop("Columns in events must be numeric.", call. = FALSE)
    }

    if (nrow(events)) {
        if (!all(is_wholenumber(events$event)))
            stop("Columns in events must be integer.", call. = FALSE)
        if (!all(is_wholenumber(events$time)))
            stop("Columns in events must be integer.", call. = FALSE)
        if (!all(is_wholenumber(events$node)))
            stop("Columns in events must be integer.", call. = FALSE)
        if (!all(is_wholenumber(events$dest)))
            stop("Columns in events must be integer.", call. = FALSE)
        if (!all(is_wholenumber(events$n)))
            stop("Columns in events must be integer.", call. = FALSE)
        if (!all(is_wholenumber(events$select)))
            stop("Columns in events must be integer.", call. = FALSE)
        if (!all(is_wholenumber(events$shift)))
            stop("Columns in events must be integer.", call. = FALSE)
    }

    events <- events[order(events$time, events$event, events$select), ]

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
              invisible(object)
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

##' Extract the events from a \code{SimInf_model} object
##'
##' Extract the scheduled events from a \code{SimInf_model} object.
##' @param model The \code{model} to extract the events from.
##' @return \code{\linkS4class{SimInf_events}} object.
##' @export
##' @examples
##' ## Create an SIR model that includes scheduled events.
##' model <- SIR(u0     = u0_SIR(),
##'              tspan  = 1:(4 * 365),
##'              events = events_SIR(),
##'              beta   = 0.16,
##'              gamma  = 0.077)
##'
##' ## Extract the scheduled events from the model and display summary
##' summary(events(model))
##'
##' ## Extract the scheduled events from the model and plot them
##' plot(events(model))
events <- function(model)
{
    check_model_argument(model)
    model@events
}

##' Extract the shift matrix from a \code{SimInf_model} object
##'
##' Utility function to extract the shift matrix \code{events@@N} from
##' a \code{SimInf_model} object, see
##' \code{\linkS4class{SimInf_events}}
##' @param model The \code{model} to extract the shift matrix
##'     \code{events@@N} from.
##' @return A mtrix.
##' @export
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Extract the shift matrix from the model
##' shift_matrix(model)
shift_matrix <- function(model)
{
    check_model_argument(model)
    model@events@N
}

##' Set the shift matrix for a \code{SimInf_model} object
##'
##' Utility function to set \code{events@@N} in a \code{SimInf_model}
##' object, see \code{\linkS4class{SimInf_events}}
##' @param model The \code{model} to set the shift matrix
##'     \code{events@@N}.
##' @param value A matrix.
##' @return \code{SimInf_model} object
##' @export
##' @importFrom methods is
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Set the shift matrix
##' shift_matrix(model) <- matrix(c(2, 1, 0), nrow = 3)
##'
##' ## Extract the shift matrix from the model
##' shift_matrix(model)
"shift_matrix<-" <- function(model, value)
{
    check_model_argument(model)

    model@events@N <- check_N(value)

    if (nrow(model@events@N) > 0 && is.null(rownames(model@events@N)))
        rownames(model@events@N) <- rownames(model@events@E)
    if (ncol(model@events@N))
        colnames(model@events@N) <- as.character(seq_len(ncol(model@events@N)))

    validObject(model)

    model
}

##' Extract the select matrix from a \code{SimInf_model} object
##'
##' Utility function to extract \code{events@@E} from a
##' \code{SimInf_model} object, see \code{\linkS4class{SimInf_events}}
##' @param model The \code{model} to extract the select matrix
##'     \code{E} from.
##' @return \code{\linkS4class{dgCMatrix}} object.
##' @export
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Extract the select matrix from the model
##' select_matrix(model)
select_matrix <- function(model)
{
    check_model_argument(model)
    model@events@E
}

##' Set the select matrix for a \code{SimInf_model} object
##'
##' Utility function to set \code{events@@E} in a \code{SimInf_model}
##' object, see \code{\linkS4class{SimInf_events}}
##' @param model The \code{model} to set the select matrix for.
##' @param value A matrix.
##' @export
##' @importFrom methods as
##' @importFrom methods is
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Set the select matrix
##' select_matrix(model) <- matrix(c(1, 0, 0, 1, 1, 1, 0, 0, 1), nrow = 3)
##'
##' ## Extract the select matrix from the model
##' select_matrix(model)
"select_matrix<-" <- function(model, value)
{
    check_model_argument(model)

    if (!is(value, "dgCMatrix"))
        value <- as(value, "dgCMatrix")

    if (!identical(Nc(model), dim(value)[1])) {
        stop("'value' must have one row for each compartment in the model.",
             call. = FALSE)
    }

    dimnames(value) <- list(rownames(model@events@E),
                            as.character(seq_len(dim(value)[2])))
    model@events@E <- value

    validObject(model)

    model
}
