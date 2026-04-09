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

degree <- function(model, a, b) {
    check_model_argument(model)

    ## Default degree is 0.
    d <- integer(n_nodes(model))

    ## Determine degree from data.
    i <- which(model@events@event == 3L)
    if (length(i) > 0) {
        dd <- tapply(methods::slot(model@events, a)[i],
                     methods::slot(model@events, b)[i],
                     function(x) {
                         length(unique(x))
                     })
        d[as.integer(dimnames(dd)[[1]])] <- dd
    }

    d
}

##' Determine in-degree for each node in a model
##'
##' Calculate the in-degree of each node based on \strong{external
##' transfer} events (\code{"extTrans"}) in the model's event
##' schedule.
##'
##' The in-degree is defined as the number of \strong{unique source
##' nodes} that have sent individuals to the target node at least
##' once.  This metric measures the connectivity of the network,
##' indicating how many different neighbors directly supply
##' individuals to a specific node.
##'
##' @param model A \code{SimInf_model} object containing the event
##'     schedule.
##' @return An integer vector where each element corresponds to a
##'     node, containing the count of unique source nodes sending
##'     individuals to it.
##' @include SimInf_model.R
##' @include check_arguments.R
##' @export
##' @examples
##' ## Create an 'SIR' model with example data.
##' model <- SIR(
##'   u0 = u0_SIR(),
##'   tspan = 1:1460,
##'   events = events_SIR(),
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Calculate the in-degree for each node.
##' deg <- indegree(model)
##'
##' ## View the in-degree for the first 6 nodes.
##' head(deg)
##'
##' ## Plot the distribution of in-degrees across all nodes.  This
##' ## shows how many source nodes typically send individuals to a given
##' ## node.
##' hist(
##'   deg,
##'   main = "Distribution of In-Degree (Unique Sources)",
##'   xlab = "Number of Unique Source Nodes"
##' )
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
