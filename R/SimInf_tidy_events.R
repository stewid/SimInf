## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2022 -- 2023 Ivana Rodriguez Ewerlöf
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

##' Class \code{"SimInf_tidy_events"}
##'
##' @slot id an integer or character identifier of the individual.
##' @slot event four event types are supported: \emph{exit},
##'     \emph{enter}, \emph{internal transfer}, and \emph{external
##'     transfer}.  When assigning the events, they can either be
##'     coded as a numerical value or a character string: \emph{exit;}
##'     \code{0} or \code{'exit'}, \emph{enter;} \code{1} or
##'     \code{'enter'}, \emph{internal transfer;} \code{2} or
##'     \code{'intTrans'}, and \emph{external transfer;} \code{3} or
##'     \code{'extTrans'}.
##' @slot time an integer, character, or date (of class \code{Date})
##'     for when the event occured. If it's a character it must be
##'     able to coerce to \code{Date}.
##' @slot node an integer or character identifier of the source node.
##' @slot dest an integer or character identifier of the destination
##'     node.
##' @export
setClass(
    "SimInf_tidy_events",
    slots = c(id    = "ANY",
              event = "integer",
              time  = "integer",
              node  = "ANY",
              dest  = "ANY"
    )
)

##' Check if a SimInf_tidy_events object is valid
##'
##' @param object The SimInf_tidy_events object to check.
##' @noRd
valid_SimInf_tidy_events <- function(object) {
    TRUE
}

## Assign the function as the validity method for the class.
setValidity("SimInf_tidy_events", valid_SimInf_tidy_events)

setAs(
    from = "SimInf_tidy_events",
    to = "data.frame",
    def = function(from) {
        events <- data.frame(id    = from@id,
                             event = from@event,
                             time  = from@time,
                             node  = from@node,
                             dest  = from@dest)

        if (!is.null(attr(events$event, "origin"))) {
            event_names <- c("exit", "enter", "intTrans", "extTrans")
            events$event <- event_names[events$event + 1]
            attr(events$event, "origin") <- NULL
        }

        if (!is.null(attr(events$time, "origin"))) {
            events$time <- as.Date(events$time,
                                   origin = attr(events$time, "origin"))
            attr(events$time, "origin") <- NULL
        }

        events
    }
)

##' Coerce to data frame
##'
##' @method as.data.frame SimInf_tidy_events
##'
##' @inheritParams base::as.data.frame
##' @export
as.data.frame.SimInf_tidy_events <- function(x, ...) {
    methods::as(x, "data.frame")
}

##' Print summary of a \code{SimInf_tidy_events} object
##'
##' @param object The \code{SimInf_tidy_events} object.
##' @return \code{invisible(object)}.
##' @export
setMethod(
    "show",
    signature(object = "SimInf_tidy_events"),
    function(object) {
        cat(sprintf(
            "Number of events: %i\n",
            length(object@id)
        ))

        invisible(object)
    }
)

##' Detailed summary of a \code{SimInf_tidy_events} object
##'
##' @param object The \code{SimInf_tidy_events} object
##' @param ... Additional arguments affecting the summary produced.
##' @return None (invisible 'NULL').
##' @export
setMethod(
    "summary",
    signature(object = "SimInf_tidy_events"),
    function(object, ...) {
        show(object)

        invisible(NULL)
    }
)

check_raw_events_id <- function(id) {
    msg <- "'id' must be an integer or character vector with non-NA values."

    if (anyNA(id))
        stop(msg, call. = FALSE)

    if (is.numeric(id)) {
        if (!all(is_wholenumber(id)))
            stop(msg, call. = FALSE)
        return(as.integer(id))
    }

    if (is.character(id))
        return(as.integer(as.factor(id)))

    if (is.factor(id))
        return(as.integer(id))

    stop(msg, call. = FALSE)
}

check_raw_events_event <- function(event) {
    msg <- "'event' must be an integer or character vector with non-NA values."

    if (anyNA(event))
        stop(msg, call. = FALSE)

    if (is.numeric(event)) {
        if (!all(is_wholenumber(event)))
            stop(msg, call. = FALSE)
        i <- as.integer(event)
        if (!all(i %in% c(0L, 1L, 3L)))
            stop(msg, call. = FALSE)
        return(i)
    }

    if (is.character(event) || is.factor(event)) {
        if (!all(event %in% c("enter", "exit", "extTrans"))) {
            stop("'event' type must be 'enter', 'exit', or 'extTrans'.",
                 call. = FALSE)
        }

        ## Find indices to 'enter', 'internal transfer' and 'external
        ## transfer' events.
        i <- rep(0L, length(event))
        i[which(event == "enter")] <- 1L
        i[which(event == "extTrans")] <- 3L
        attr(i, "origin") <- "character"
        return(i)
    }

    stop(msg, call. = FALSE)
}

check_raw_events_time <- function(time) {
    msg <- "'time' must be an integer or character vector with non-NA values."

    if (is.numeric(time)) {
        if (anyNA(time))
            stop(msg, call. = FALSE)
        if (!all(is_wholenumber(time)))
            stop(msg, call. = FALSE)
        return(as.integer(time))
    }

    if (is.character(time) || is.factor(time))
        time <- as.Date(time)

    if (anyNA(time))
        stop(msg, call. = FALSE)

    if (methods::is(time, "Date")) {
        time <- as.integer(julian(time, origin = as.Date("1970-01-01")))
        attr(time, "origin") <- "1970-01-01"
        return(time)
    }

    stop(msg, call. = FALSE)
}

check_raw_events_nodes <- function(event, node, dest) {
    if (any(anyNA(node), anyNA(dest[event == 3L]))) {
        stop("'node' or 'dest' contain NA values.",
             call. = FALSE)
    }

    if (all(is.numeric(node), is.numeric(dest))) {
        if (any(!all(is_wholenumber(node)),
                !all(is_wholenumber(dest[event == 3L])))) {
            stop("'node' and 'dest' must both be integer or character.",
                 call. = FALSE)
        }

        node <- as.integer(node)
        dest <- as.integer(dest)

        return(list(node = node, dest = dest))
    }

    if (all(is.character(node), is.character(dest))) {
        nodes <- as.factor(unique(c(node, dest)))
        node <- as.integer(factor(node, levels = levels(nodes)))
        dest <- as.integer(factor(dest, levels = levels(nodes)))
        return(list(node = node, dest = dest))
    }

    stop("'node' and 'dest' must both be integer or character.",
         call. = FALSE)
}

##' Tidy events
##'
##' In many countries, individual-based livestock data are collected
##' to enable contact tracing during disease outbreaks. However, the
##' livestock databases are not always structured in such a way that
##' relevant information for disease spread simulations is easily
##' retrieved. The aim of this function is to facilitate cleaning
##' livestock event data and prepare it for usage in SimInf.
##'
##' The argument \code{events} in \code{tidy_events} must be a
##' \code{data.frame} with the following columns:
##' * **id:** an integer or character identifier of the individual.
##' * **event:** four event types are supported: \emph{exit},
##'     \emph{enter}, \emph{internal transfer}, and \emph{external
##'     transfer}.  When assigning the events, they can either be
##'     coded as a numerical value or a character string: \emph{exit;}
##'     \code{0} or \code{'exit'}, \emph{enter;} \code{1} or
##'     \code{'enter'}, \emph{internal transfer;} \code{2} or
##'     \code{'intTrans'}, and \emph{external transfer;} \code{3} or
##'     \code{'extTrans'}.
##' * **time:** an integer, character, or date (of class \code{Date})
##'     for when the event occured. If it's a character it must be
##'     able to coerce to \code{Date}.
##' * **node:** an integer or character identifier of the source node.
##' * **dest:** an integer or character identifier of the destination
##'     node for movement events, and 'dest' will be replaced with
##'     \code{NA} for non-movement events.
##' @param events a \code{data.frame} with the columns `id`, `event`,
##'     `time`, `node`, and `dest` to define the events, see
##'     `details`.
##' @return \linkS4class{SimInf_tidy_events}
##' @export
##' @md
tidy_events <- function(events) {
    columns <- c("id", "event", "time", "node", "dest")
    if (!is.data.frame(events))
        events <- as.data.frame(events)
    if (!all(columns %in% names(events)))
        stop("Missing columns in 'events'.", call. = FALSE)
    events <- events[, columns, drop = FALSE]

    id <- check_raw_events_id(events$id)
    event <- check_raw_events_event(events$event)
    time <- check_raw_events_time(events$time)
    events$dest[event != 3L] <- NA
    nodes <- check_raw_events_nodes(event, events$node, events$dest)

    keep <- .Call(SimInf_clean_raw_events,
                  id,
                  event,
                  time,
                  nodes$node,
                  nodes$dest)

    origin <- attr(time, "origin")
    time <- time[keep]
    if (!is.null(origin))
        attr(time, "origin") <- origin

    origin <- attr(event, "origin")
    event <- event[keep]
    if (!is.null(origin))
        attr(event, "origin") <- origin

    methods::new("SimInf_tidy_events",
                 id    = events$id[keep],
                 event = event,
                 time  = time,
                 node  = events$node[keep],
                 dest  = events$dest[keep])
}

## Check for a valid 'at' parameter
tidy_events_at <- function(events, at) {
    if (is.null(at))
        return(min(events@time))

    if (is.numeric(at)) {
        if (!all(is_wholenumber(at)))
            stop("'at' must be an integer or date.", call. = FALSE)
        return(as.integer(at))
    }

    stop("Not implemented.")
}

##' Extract individuals from tidy events
##'
##' Lookup individuals, in which node they are located and their age
##' at a specified time-point \code{i}.
##' @rdname SimInf_tidy_events-index-methods
##' @param x a tidy events \code{object}.
##' @param i FIXME.
##' @return a \code{data.frame} with the columns \code{id},
##'     \code{node}, and \code{age}.
##' @export
setMethod(
    "[",
    signature(x = "SimInf_tidy_events", i = "ANY", j = "ANY"),
    function(x, i) {
        if (missing(i))
            i <- NULL

        ## Check that all individuals have an enter event.
        if (length(setdiff(x@id, x@id[x@event == 1L])))
            stop("All individuals must have an 'enter' event.", call. = FALSE)

        ## Keep events for 'u0' that are <= 'at'. Drop individuals
        ## that exit before 'at'.
        k <- which(x@time <= tidy_events_at(x, i))
        drop <- unique(x@id[k[x@event[k] == 0L]])
        k <- k[!(x@id[k] %in% drop)]

        ## Keep last event for each individual.
        l <- as.integer(tapply(k, x@id[k], max))

        ## If it's a movement, swap node and dest to have the current
        ## node location of the individual.
        node <- ifelse(x@event[l] == 3L, x@dest[l], x@node[l])

        ## Determine age in days of each individual by subtracting the
        ## current time with the first event for each individual.
        age <- x@time[l] - x@time[as.integer(tapply(k, x@id[k], min))]

        data.frame(id = x@id[l], node = node, age = age)
    }
)

##' Display the distribution of raw events over time
##'
##' @param x The raw events data to plot.
##' @param frame.plot Draw a frame around each plot. Default is FALSE.
##' @param ... Additional arguments affecting the plot
##' @aliases plot,SimInf_tidy_events-method
##' @export
setMethod(
    "plot",
    signature(x = "SimInf_tidy_events"),
    function(x, frame.plot = FALSE, ...) {
        savepar <- graphics::par(mfrow = c(2, 2),
                                 oma = c(1, 1, 2, 0),
                                 mar = c(4, 3, 1, 1))
        on.exit(graphics::par(savepar))

        n <- rep(1L, length(x@id))
        yy <- stats::xtabs(n ~ event + time,
                           cbind(event = x@event,
                                 time = x@time,
                                 n = n))

        xx <- as.integer(colnames(yy))
        if (!is.null(attr(x@time, "origin")))
            xx <- as.Date(xx, origin = attr(x@time, "origin"))

        ## Plot events
        plot_SimInf_events(xx, yy, "Exit", frame.plot, ...)
        plot_SimInf_events(xx, yy, "Enter", frame.plot, ...)
        plot_SimInf_events(xx, yy, "Internal transfer", frame.plot, ...)
        plot_SimInf_events(xx, yy, "External transfer", frame.plot, ...)
    }
)
