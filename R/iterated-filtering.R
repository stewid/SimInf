## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2025 Stefan Widgren
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

##' Class \code{"SimInf_if2"}
##'
##' @slot model The \code{SimInf_model} object to estimate parameters
##'     in.
##' @template priors-slot
##' @include SimInf_model.R
##' @export
setClass(
    "SimInf_if2",
    slots = c(model = "SimInf_model",
              priors = "data.frame")
)

##' Brief summary of a \code{SimInf_if2} object
##'
##' @param object The \code{SimInf_if2} object.
##' @return \code{invisible(object)}.
##' @export
setMethod(
    "show",
    signature(object = "SimInf_if2"),
    function(object) {
        cat("Iterated filtering\n")
        cat("------------------\n")

        invisible(object)
    }
)

##' Detailed summary of a \code{SimInf_if2} object
##'
##' @param object The \code{SimInf_if2} object
##' @param ... Not used.
##' @return None (invisible 'NULL').
##' @export
setMethod(
    "summary",
    signature(object = "SimInf_if2"),
    function(object, ...) {
        cat("Iterated filtering\n")
        cat("------------------\n")

        invisible(NULL)
    }
)
