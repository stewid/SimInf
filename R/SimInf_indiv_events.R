## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2022 -- 2023 Ivana Rodriguez Ewerlöf
## Copyright (C) 2015 -- 2024 Stefan Widgren
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

##' Class \code{"SimInf_indiv_events"}
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
    "SimInf_indiv_events",
    slots = c(id    = "ANY",
              event = "integer",
              time  = "integer",
              node  = "ANY",
              dest  = "ANY"
    )
)

##' Check if a SimInf_indiv_events object is valid
##'
##' @param object The SimInf_indiv_events object to check.
##' @noRd
valid_SimInf_indiv_events <- function(object) {
    TRUE
}

## Assign the function as the validity method for the class.
setValidity("SimInf_indiv_events", valid_SimInf_indiv_events)

setAs(
    from = "SimInf_indiv_events",
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
##' @method as.data.frame SimInf_indiv_events
##'
##' @inheritParams base::as.data.frame
##' @export
as.data.frame.SimInf_indiv_events <- function(x, ...) {
    methods::as(x, "data.frame")
}

##' Print summary of a \code{SimInf_indiv_events} object
##'
##' @param object The \code{SimInf_indiv_events} object.
##' @return \code{invisible(object)}.
##' @export
setMethod(
    "show",
    signature(object = "SimInf_indiv_events"),
    function(object) {
        cat(sprintf(
            paste0("Number of individuals: %i\n",
                   "Number of events: %i\n"),
            length(unique(object@id)),
            length(object@id)
        ))

        invisible(object)
    }
)

##' Detailed summary of a \code{SimInf_indiv_events} object
##'
##' @param object The \code{SimInf_indiv_events} object
##' @param ... Additional arguments affecting the summary produced.
##' @return None (invisible 'NULL').
##' @export
setMethod(
    "summary",
    signature(object = "SimInf_indiv_events"),
    function(object, ...) {
        show(object)

        for (i in seq_len(4)) {
            switch(i,
                   cat(" - Exit: "),
                   cat(" - Enter: "),
                   cat(" - Internal transfer: "),
                   cat(" - External transfer: "))

            j <- which(object@event == (i - 1L))
            cat(sprintf("%i\n", length(j)))
        }

        invisible(NULL)
    }
)

check_indiv_events_id <- function(id) {
    msg <- "'id' must be an integer or character vector with non-NA values."

    if (anyNA(id))
        stop(msg, call. = FALSE)

    if (is.numeric(id)) {
        if (!all(is_wholenumber(id)))
            stop(msg, call. = FALSE)
        return(as.integer(id))
    }

    if (is.character(id))
        id <- as.factor(id)

    if (is.factor(id))
        return(as.integer(id))

    stop(msg, call. = FALSE)
}

check_indiv_events_event <- function(event) {
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

check_indiv_events_time <- function(time) {
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

check_indiv_events_nodes <- function(event, node, dest) {
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

##' Individual events
##'
##' In many countries, individual-based livestock data are collected
##' to enable contact tracing during disease outbreaks. However, the
##' livestock databases are not always structured in such a way that
##' relevant information for disease spread simulations is easily
##' retrieved. The aim of this function is to facilitate cleaning
##' livestock event data and prepare it for usage in SimInf.
##'
##' The argument \code{events} in \code{individual_events} must be a
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
##' @return \linkS4class{SimInf_indiv_events}
##' @export
##' @md
individual_events <- function(events) {
    columns <- c("id", "event", "time", "node", "dest")
    if (!is.data.frame(events))
        events <- as.data.frame(events)
    if (!all(columns %in% names(events)))
        stop("Missing columns in 'events'.", call. = FALSE)
    events <- events[, columns, drop = FALSE]

    id <- check_indiv_events_id(events$id)
    event <- check_indiv_events_event(events$event)
    time <- check_indiv_events_time(events$time)
    events$dest[event != 3L] <- NA
    nodes <- check_indiv_events_nodes(event, events$node, events$dest)

    ## Ensure the events are sorted.
    i <- order(id, time, event)

    keep <- .Call(SimInf_clean_indiv_events,
                  id[i],
                  event[i],
                  time[i],
                  nodes$node[i],
                  nodes$dest[i])

    origin <- attr(time, "origin")
    time <- time[i][keep]
    if (!is.null(origin))
        attr(time, "origin") <- origin

    origin <- attr(event, "origin")
    event <- event[i][keep]
    if (!is.null(origin))
        attr(event, "origin") <- origin

    methods::new("SimInf_indiv_events",
                 id    = events$id[i][keep],
                 event = event,
                 time  = time,
                 node  = events$node[i][keep],
                 dest  = events$dest[i][keep])
}

## Check for a valid 'time' parameter
indiv_events_time <- function(events, time) {
    msg <- "'time' must be an integer or date."

    if (is.null(time))
        return(min(events@time))

    if (is.numeric(time)) {
        if (any(!all(is_wholenumber(time)), length(time) != 1L, anyNA(time)))
            stop(msg, call. = FALSE)
        return(as.integer(time))
    }

    if (is.character(time) || is.factor(time))
        time <- as.Date(time)

    if (any(anyNA(time), length(time) != 1L))
        stop(msg, call. = FALSE)

    if (methods::is(time, "Date")) {
        origin <- attr(events@time, "origin")
        if (is.null(origin))
            stop("'time' must be an integer.", call. = FALSE)
        origin <- as.Date(origin)
        return(as.integer(julian(time, origin = origin)))
    }

    stop(msg, call. = FALSE)
}

##' Extract individuals from \code{SimInf_indiv_events}
##'
##' Lookup individuals, in which node they are located and their age
##' at a specified time-point.
##' @param x an individual events object of class
##'     \code{SimInf_indiv_events}.
##' @param time the time-point for the lookup of individuals. Default
##'     is \code{NULL} which means to extract the individuals at the
##'     minimum time-point in the events object \code{x}.
##' @return a \code{data.frame} with the columns \code{id},
##'     \code{node}, and \code{age}.
##' @export
setGeneric(
    "get_individuals",
    signature = "x",
    function(x, time = NULL) {
        standardGeneric("get_individuals")
    }
)

##' @rdname get_individuals
##' @export
setMethod(
    "get_individuals",
    signature(x = "SimInf_indiv_events"),
    function(x, time = NULL) {
        ## Check that all individuals have an enter event.
        if (length(setdiff(x@id, x@id[x@event == 1L])))
            stop("All individuals must have an 'enter' event.", call. = FALSE)

        ## Keep events that occur <= 'time'. Drop individuals that
        ## exit before 'time'.
        time <- indiv_events_time(x, time)
        k <- which(x@time <= time)
        drop <- unique(x@id[k[x@event[k] == 0L]])
        k <- k[!(x@id[k] %in% drop)]

        ## Keep last event for each individual.
        l <- as.integer(tapply(k, x@id[k], max))

        ## If it's a movement, swap node and dest to have the current
        ## node location of the individual.
        node <- ifelse(x@event[l] == 3L, x@dest[l], x@node[l])

        ## Determine age in days of each individual by subtracting the
        ## current time with the first event for each individual.
        age <- time - x@time[as.integer(tapply(k, x@id[k], min))]

        data.frame(id = x@id[l], node = node, age = age)
    }
)

##' Display the distribution of individual events over time
##'
##' @param x The individual events data to plot.
##' @template plot-frame-param
##' @param ... Other graphical parameters that are passed on to the
##'     plot function.
##' @aliases plot,SimInf_indiv_events-method
##' @export
setMethod(
    "plot",
    signature(x = "SimInf_indiv_events"),
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

##' @rdname events
##' @param time Only used when object is of class
##'     \code{SimInf_indiv_events} object. All events that occur after
##'     \sQuote{time} are included. Default is \code{NULL} which means
##'     to extract the events after the minimum time-point in the
##'     \code{SimInf_indiv_events} object.
##' @param target Only used when object is of class
##'     \code{SimInf_indiv_events} object. The SimInf model ('SEIR',
##'     'SIR', 'SIS', 'SISe3', 'SISe3_sp', 'SISe', or 'SISe_sp') to
##'     target the events and u0 for. The default, \code{NULL},
##'     creates events but they might have to be post-processed to fit
##'     the specific use case.
##' @param age Only used when object is of class
##'     \code{SimInf_indiv_events} object. Integer vector with break
##'     points in days for the ageing events.
##' @export
setMethod(
    "events",
    signature(object = "SimInf_indiv_events"),
    function(object, time = NULL, target = NULL, age = NULL) {
        ## Check for valid target model.
        if (!is.null(target)) {
            target <- match.arg(target, c("SEIR", "SIR", "SIS",
                                          "SISe3", "SISe3_sp", "SISe",
                                          "SISe_sp"))
        } else {
            stop("Not implemented.", call. = FALSE)
        }

        ## Keep events that occur after 'time'.
        i <- which(object@time > indiv_events_time(object, time))

        node <- object@node[i]
        dest <- object@dest[i]
        dest[is.na(dest)] <- 0L

        data.frame(event      = object@event[i],
                   time       = object@time[i],
                   node       = node,
                   dest       = dest,
                   n          = 1,
                   proportion = 0,
                   select     = 1,
                   shift      = 0)
    }
)
