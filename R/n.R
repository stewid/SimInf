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

##' Determine the number of nodes in a model
##'
##' Extract the number of nodes from a \code{SimInf_model} object.  A
##' node represents a distinct sub-population in the model, but its
##' definition is determined by the modeller. For example, a node can
##' represent a cattle herd, a pen within a herd, or even a single
##' individual, depending on the research question and the scale of
##' the study. This count is equivalent to the number of rows in the
##' initial state vector \code{u0}.
##'
##' @param model A \code{SimInf_model} object.
##' @return An integer scalar representing the total number of nodes
##'     in the model.
##' @export
##' @examples
##' ## Create an 'SIR' model with 100 nodes.
##' u0 <- data.frame(
##'   S = rep(99, 100),
##'   I = rep(1, 100),
##'   R = rep(0, 100)
##' )
##'
##' model <- SIR(
##'   u0 = u0,
##'   tspan = 1:10,
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Get the number of nodes.
##' n_nodes(model)
## nolint start: brace_linter
setGeneric(
    "n_nodes",
    signature = "model",
    function(model)
        standardGeneric("n_nodes")
)
## nolint end

##' @rdname n_nodes
##' @include SimInf_model.R
##' @export
setMethod(
    "n_nodes",
    signature(model = "SimInf_model"),
    function(model) {
        as.integer(floor(dim(model@u0)[2] / n_replicates(model)))
    }
)

##' @rdname n_nodes
##' @include pfilter.R
##' @export
setMethod(
    "n_nodes",
    signature(model = "SimInf_pfilter"),
    function(model) {
        n_nodes(model@model)
    }
)

##' @rdname n_nodes
##' @include pmcmc.R
##' @export
setMethod(
    "n_nodes",
    signature(model = "SimInf_pmcmc"),
    function(model) {
        n_nodes(model@model)
    }
)

##' Determine the number of replicates in a model
##'
##' Extract the number of replicates from a \code{SimInf_model}
##' object.  Replicates are independent copies of the model state used
##' by the \strong{particle filter} (\code{\link{pfilter}}).  This
##' value is set \strong{internally} by \code{pfilter} and is not
##' specified when creating a standard model. If the model has not
##' been processed by a filter, the value will be 1.
##'
##' @param model A \code{SimInf_model} object.
##' @return An integer scalar representing the number of replicates
##'     in the model.
##' @export
##' @examples
##' ## Create a standard 'SIR' model.
##' u0 <- data.frame(
##'   S = rep(99, 100),
##'   I = rep(1, 100),
##'   R = rep(0, 100)
##' )
##'
##' model <- SIR(
##'   u0 = u0,
##'   tspan = 1:10,
##'   beta = 0.16,
##'   gamma = 0.077
##' )
##'
##' ## Get the number of replicates (default is 1 for standard
##' ## models).
##' n_replicates(model)
## nolint start: brace_linter
setGeneric(
    "n_replicates",
    signature = "model",
    function(model)
        standardGeneric("n_replicates")
)
## nolint end

##' @rdname n_replicates
##' @include SimInf_model.R
##' @export
setMethod(
    "n_replicates",
    signature(model = "SimInf_model"),
    function(model) {
        model@replicates
    }
)

##' Determine the number of compartments in a model
##'
##' @param model the \code{model} object to extract the number of
##'     compartments from.
##' @return the number of compartments in the model.
##' @export
##' @examples
##' ## Create an 'SIR' model with 100 nodes, with 99 susceptible,
##' ## 1 infected and 0 recovered in each node.
##' u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Display the number of compartments in the model.
##' n_compartments(model)
## nolint start: brace_linter
setGeneric(
    "n_compartments",
    signature = "model",
    function(model)
        standardGeneric("n_compartments")
)
## nolint end

##' @rdname n_compartments
##' @include SimInf_model.R
##' @export
setMethod(
    "n_compartments",
    signature(model = "SimInf_model"),
    function(model) {
        dim(model@S)[1]
    }
)

## Number of compartments
Nc <- function(model) {
    n_compartments(model)
}

##' Determine the number of transitions in a model
##'
##' @param model the \code{model} object to extract the number of
##'     transitions from.
##' @return the number of transitions in the model.
##' @noRd
##' @examples
##' ## Create an 'SIR' model with 100 nodes, with 99 susceptible,
##' ## 1 infected and 0 recovered in each node.
##' u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Display the number of transitions in the model.
##' n_transitions(model)
## nolint start: brace_linter
setGeneric(
    "n_transitions",
    signature = "model",
    function(model)
        standardGeneric("n_transitions")
)
## nolint end

##' @rdname n_transitions
##' @include SimInf_model.R
##' @noRd
setMethod(
    "n_transitions",
    signature(model = "SimInf_model"),
    function(model) {
        dim(model@G)[1]
    }
)

## Number of continuous state variables
Nd <- function(model) {
    check_model_argument(model)
    dim(model@v0)[1]
}
