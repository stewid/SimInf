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

##' @importFrom methods slot
##' @noRd
degree <- function(model, a, b) {
    check_model_argument(model)

    ## Default degree is 0.
    d <- integer(Nn(model))

    ## Determine degree from data.
    i <- which(model@events@event == 3L)
    if (length(i) > 0) {
        dd <- tapply(slot(model@events, a)[i],
                     slot(model@events, b)[i],
                     function(x) {
                         length(unique(x))
                     })
        d[as.integer(dimnames(dd)[[1]])] <- dd
    }

    d
}

##' Determine in-degree for each node in a model
##'
##' The number of nodes with inward \emph{external transfer} events to
##' each node.
##' @param model determine in-degree for each node in the model.
##' @return vector with in-degree for each node.
##' @include SimInf_model.R
##' @include check_arguments.R
##' @export
##' @examples
##' ## Create an 'SIR' model with 1600 nodes and initialize
##' ## it with example data.
##' model <- SIR(u0 = u0_SIR(), tspan = 1:1460, events = events_SIR(),
##'              beta   = 0.16, gamma  = 0.077)
##'
##' ## Display indegree for each node in the model.
##' plot(indegree(model))
indegree <- function(model) {
    degree(model, "node", "dest")
}

##' Determine out-degree for each node in a model
##'
##' The number nodes that are connected with \emph{external transfer}
##' events from each node.
##' @param model determine out-degree for each node in the model.
##' @return vector with out-degree for each node.
##' @include SimInf_model.R
##' @include check_arguments.R
##' @export
##' @examples
##' ## Create an 'SIR' model with 1600 nodes and initialize
##' ## it with example data.
##' model <- SIR(u0 = u0_SIR(), tspan = 1:1460, events = events_SIR(),
##'              beta   = 0.16, gamma  = 0.077)
##'
##' ## Display outdegree for each node in the model.
##' plot(outdegree(model))
outdegree <- function(model) {
    degree(model, "dest", "node")
}
