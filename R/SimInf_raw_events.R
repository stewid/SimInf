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

##' Class \code{"SimInf_raw_events"}
##'
##' @slot id an integer or character identifier of the source holding
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
##' @slot keep logical vector with the events that are kept after the
##'     cleaning.
##' @export
setClass(
    "SimInf_raw_events",
    slots = c(id    = "ANY",
              event = "ANY",
              time  = "ANY",
              node  = "ANY",
              dest  = "ANY",
              keep  = "logical"
    )
)

##' Check if a SimInf_raw_events object is valid
##'
##' @param object The SimInf_raw_events object to check.
##' @noRd
valid_SimInf_raw_events_object <- function(object) {
    TRUE
}

## Assign the function as the validity method for the class.
setValidity("SimInf_raw_events", valid_SimInf_raw_events_object)

setAs(
    from = "SimInf_raw_events",
    to = "data.frame",
    def = function(from) {
        data.frame(id    = from@id[from@keep],
                   event = from@event[from@keep],
                   time  = from@time[from@keep],
                   node  = from@node[from@keep],
                   dest  = from@dest[from@keep]
        )
    }
)

##' Coerce to data frame
##'
##' @method as.data.frame SimInf_raw_events
##'
##' @inheritParams base::as.data.frame
##' @export
as.data.frame.SimInf_raw_events <- function(x, ...) {
    methods::as(x, "data.frame")
}

##' Print summary of a \code{SimInf_raw_events} object
##'
##' @param object The \code{SimInf_raw_events} object.
##' @return \code{invisible(object)}.
##' @export
setMethod(
    "show",
    signature(object = "SimInf_raw_events"),
    function(object) {
        cat(sprintf(
            "Number of events: %i\n",
            length(object@id)
        ))
        cat(sprintf(
            "Proportion of events to keep: %f\n",
            mean(object@keep)
        ))

        invisible(object)
    }
)

##' Detailed summary of a \code{SimInf_raw_events} object
##'
##' @param object The \code{SimInf_raw_events} object
##' @param ... Additional arguments affecting the summary produced.
##' @return None (invisible 'NULL').
##' @export
setMethod(
    "summary",
    signature(object = "SimInf_raw_events"),
    function(object, ...) {
        show(object)

        invisible(NULL)
    }
)

##' Raw events
##'
##' @param events FIXME
##' @return FIXME
##' @export
raw_events <- function(events) {
    events <- check_raw_events(events)

    keep <- .Call(SimInf_clean_raw_events,
                  events$id,
                  events$event,
                  events$time,
                  events$node,
                  events$dest)

    methods::new("SimInf_raw_events",
                 id    = events$id,
                 event = events$event,
                 time  = events$time,
                 node  = events$node,
                 dest  = events$dest,
                 keep  = keep)
}

##' Display the distribution of raw events over time
##'
##' @param x The raw events data to plot.
##' @param frame.plot Draw a frame around each plot. Default is FALSE.
##' @param ... Additional arguments affecting the plot
##' @aliases plot,SimInf_raw_events-method
##' @export
setMethod(
    "plot",
    signature(x = "SimInf_raw_events"),
    function(x, frame.plot = FALSE, ...) {
        savepar <- graphics::par(mfrow = c(2, 2),
                                 oma = c(1, 1, 2, 0),
                                 mar = c(4, 3, 1, 1))
        on.exit(graphics::par(savepar))

        n <- rep(1L, length(x@keep))
        all_y <- stats::xtabs(n ~ event + time,
                           cbind(event = x@event,
                                 time = x@time,
                                 n = n))
        all_x <- as.integer(colnames(all_y))
        if (!is.null(attr(x@time, "origin")))
            all_x <- as.Date(all_x, origin = attr(x@time, "origin"))

        keep_y <- stats::xtabs(n ~ event + time,
                           cbind(event = x@event[x@keep],
                                 time = x@time[x@keep],
                                 n = n[x@keep]))
        keep_x <- as.integer(colnames(keep_y))
        if (!is.null(attr(x@time, "origin")))
            keep_x <- as.Date(keep_x, origin = attr(x@time, "origin"))

        ## Plot events
        plot_SimInf_events(all_x, all_y, "Exit", frame.plot,
                           col = "gray70", ...)
        if (length(keep_x) && "0" %in% rownames(keep_y))
            lines(keep_x, keep_y["0", ])

        plot_SimInf_events(all_x, all_y, "Enter", frame.plot,
                           col = "gray70", ...)
        if (length(keep_x) && "1" %in% rownames(keep_y))
            lines(keep_x, keep_y["1", ])

        plot_SimInf_events(all_x, all_y, "Internal transfer",
                           frame.plot, col = "gray70", ...)
        if (length(keep_x) && "2" %in% rownames(keep_y))
            lines(keep_x, keep_y["2", ])

        plot_SimInf_events(all_x, all_y, "External transfer",
                           frame.plot, col = "gray70", ...)
        if (length(keep_x) && "3" %in% rownames(keep_y))
            lines(keep_x, keep_y["3", ])
    }
)
