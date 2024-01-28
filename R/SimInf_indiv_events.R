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

inject_ageing_events <- function(events, age) {
    ## Create ageing events for individuals. First, determine the
    ## time-points for the enter events.
    ageing <- events[events[, "event"] == 1L, c("id", "time"), drop = FALSE]
    colnames(ageing) <- c("id", "enter")

    ## Then determine the time-points for the exit events.
    i <- which(events[, "event"] == 0L)
    j <- match(ageing[, "id"], events[i, "id"])
    ageing <- cbind(ageing, exit = events[i, "time"][j])

    ## Ensure all ageing events occur within the time-span of all
    ## events.
    ageing[is.na(ageing[, "exit"]), "exit"] <- max(events[, "time"]) + 1L

    ## Add a column of the age in days of the individual at
    ## each event.
    i <- match(events[, "id"], ageing[, "id"])
    events <- cbind(events, days = events[, "time"] - ageing[i, "enter"])

    ## Clear select for non-enter events.
    events[events[, "event"] != 1L, "select"] <- NA_integer_

    for (i in seq_len(length(age))[-1]) {
        ## Determine which individuals are eligble for ageing.
        j <- which(ageing[, "exit"] - ageing[, "enter"] > age[i])
        ageing <- ageing[j, , drop = FALSE]

        if (nrow(ageing) > 0) {
            events <- rbind(
                events,
                matrix(c(ageing[, "id"],
                         rep(2L, nrow(ageing)),
                         ageing[, "enter"] + age[i],
                         rep(0L, nrow(ageing)),
                         rep(0L, nrow(ageing)),
                         rep(length(age) + i - 1L, nrow(ageing)),
                         rep(i - 1L, nrow(ageing)),
                         rep(age[i], nrow(ageing))),
                       nrow = nrow(ageing),
                       ncol = 8L,
                       dimnames = list(
                           NULL,
                           c("id", "event", "time", "node", "dest",
                             "select", "shift", "days"))))

            events <- events[order(events[, "id"],
                                   events[, "time"],
                                   events[, "event"],
                                   events[, "node"],
                                   events[, "dest"]), ]

            ## Determine the node for the ageing event.
            j <- which(events[, "event"] == 2L &
                       events[, "node"] == 0L)
            if (length(j)) {
                ## If the previous event is a movement, then the node
                ## is the destination of that movement.
                events[j, "node"] <- ifelse(events[j - 1, "dest"] > 0L,
                                            events[j - 1, "dest"],
                                            events[j - 1, "node"])
            }
        }

        ## Determine if there are any exit events for this age
        ## category.
        j <- which(events[, "event"] == 0L &
                   events[, "days"] <= age[i] &
                   is.na(events[, "select"]))
        events[j, "select"] <- length(age) + i - 1L

        ## Determine if there are any movement events for this age
        ## category.
        j <- which(events[, "event"] == 3L &
                   events[, "days"] < age[i] &
                   is.na(events[, "select"]))
        events[j, "select"] <- length(age) + i - 1L
    }

    ## Determine if there are any events where select = NA because of
    ## individuals being older than the largest age break point.
    j <- which(is.na(events[, "select"]))
    events[j, "select"] <- 2L * length(age)

    events
}

events_target <- function(events, target) {
    if (!is.null(target)) {
        if (target %in% c("SIS", "SISe", "SISe_sp", "SEIR")) {
            events[events[, "event"] == 0L, "select"] <- 2L
            events[events[, "event"] == 3L, "select"] <- 2L
        } else if (target %in% c("SIR")) {
            events[events[, "event"] == 0L, "select"] <- 4L
            events[events[, "event"] == 3L, "select"] <- 4L
        }
    }

    events <- as.data.frame(events)
    events$proportion <- as.numeric(events$proportion)
    events
}

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
##' @section Transform individual events:

##'     In many countries, individual-based livestock data are
##'     collected to enable contact tracing during disease
##'     outbreaks. However, the livestock databases are not always
##'     structured in such a way that relevant information for disease
##'     spread simulations is easily retrieved. The aim of this
##'     function is to facilitate cleaning livestock event data and
##'     prepare it for usage in SimInf.
##' @export
setMethod(
    "events",
    signature(object = "SimInf_indiv_events"),
    function(object, time = NULL, target = NULL, age = NULL) {
        age <- check_age(age)
        target <- check_target(target, age)

        ## Map nodes to the one-based index in SimInf.
        all_nodes <- unique(c(object@node, object@dest))
        all_nodes <- sort(all_nodes[!is.na(all_nodes)])

        events <- matrix(c(as.integer(as.factor(object@id)),
                           as.integer(object@event),
                           as.integer(object@time),
                           match(object@node, all_nodes),
                           match(object@dest, all_nodes),
                           rep(1L, length(object@id)),
                           rep(0L, length(object@id))),
                         nrow = length(object@id),
                         ncol = 7,
                         dimnames = list(
                             NULL,
                             c("id", "event", "time", "node", "dest",
                               "select", "shift")))

        ## Ensure all 'dest' are non-NA.
        events[is.na(events[, "dest"]), "dest"] <- 0L

        ## Check that all individuals have an enter event.
        if (length(setdiff(events[, "id"],
                           events[events[, "event"] == 1L, "id"]))) {
            stop("All individuals must have an 'enter' event.", call. = FALSE)
        }

        if (length(age) > 1)
            events <- inject_ageing_events(events, age)

        ## Keep events that occur after 'time'.
        i <- which(events[, "time"] > indiv_events_time(object, time))
        j <- order(events[i, "time"], events[i, "event"], events[i, "node"],
                   events[i, "dest"], events[i, "select"], events[i, "shift"])
        events <- cbind(events[i[j], c("event", "time", "node", "dest", "select",
                                       "shift"), drop = FALSE],
                        n = 1L)

        if (nrow(events) > 1) {
            ## Determine duplicated events.
            i <- seq(from = 2L, to = nrow(events), by = 1L)
            duplicated <- i[events[i, "time"] == events[i - 1, "time"] &
                            events[i, "event"] == events[i - 1, "event"] &
                            events[i, "node"] == events[i - 1, "node"] &
                            events[i, "dest"] == events[i - 1, "dest"] &
                            events[i, "select"] == events[i - 1, "select"] &
                            events[i, "shift"] == events[i - 1, "shift"]]

            if (length(duplicated) > 0) {
                ## Determine the index to the first non-duplicated row
                ## before each sequence of duplicates.
                d <- diff(duplicated)
                i <- duplicated[c(0L, d) != 1L] - 1L

                ## Count the number of duplicates in each sequence.
                n <- cumsum(c(1L, d) > 1L)
                n <- as.integer(tapply(n, n, length))

                ## Update events and remove duplicates.
                events[i, "n"] <- events[i, "n"] + n
                events <- events[-duplicated, ]
            }
        }

        events <- cbind(events, proportion = 0L)
        rownames(events) <- NULL

        events <- events_target(events[, c("event", "time", "node",
                                           "dest", "n", "proportion",
                                           "select", "shift"),
                                       drop = FALSE],
                                target)

        if (!is.null(attr(object@event, "origin"))) {
            event_names <- c("exit", "enter", "intTrans", "extTrans")
            events$event <- event_names[events$event + 1]
        }

        if (!is.null(attr(object@time, "origin"))) {
            events$time <- as.Date(events$time,
                                   origin = attr(object@time, "origin"))
            events$time <- as.character(events$time)
        } else {
            events$time <- as.numeric(events$time)
        }

        events
    }
)
