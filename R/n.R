## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2020 Stefan Widgren
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
##' @param model the \code{model} object to extract the number of
##'     nodes from.
##' @return the number of nodes in the model.
##' @export
##' @examples
##' ## Create an 'SIR' model with 100 nodes, with 99 susceptible,
##' ## 1 infected and 0 recovered in each node.
##' u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Display the number of nodes in the model.
##' n_nodes(model)
setGeneric(
    "n_nodes",
    signature = "model",
    function(model) {
        standardGeneric("n_nodes")
    }
)

##' @rdname n_nodes
##' @include SimInf_model.R
##' @export
setMethod(
    "n_nodes",
    signature(model = "SimInf_model"),
    function(model) {
        dim(model@u0)[2]
    }
)

##' Determine the number of nodes in a model
##'
##' Extract number of nodes in a model.
##' @param model the \code{model} object to extract the number of
##'     nodes from.
##' @return the number of nodes in the model.
##' @export
##' @examples
##' ## Create an 'SIR' model with 100 nodes, with 99 susceptible,
##' ## 1 infected and 0 recovered in each node.
##' u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Display the number of nodes in the model.
##' Nn(model)
Nn <- function(model) {
    .Deprecated(new = "Use 'n_nodes' instead.",
                msg = paste("'Nn' will be removed in a future version,",
                            "use 'n_nodes' instead."))
    n_nodes(model)
}

## Number of compartments
Nc <- function(model) {
    check_model_argument(model)
    dim(model@S)[1]
}

## Number of transitions
Nt <- function(model) {
    check_model_argument(model)
    dim(model@G)[1]
}

## Number of continuous state variables
Nd <- function(model) {
    check_model_argument(model)
    dim(model@v0)[1]
}
