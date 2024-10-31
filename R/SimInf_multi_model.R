## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
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

##' Class \code{"SimInf_multi_model"}
##'
##' @slot model The \code{SimInf_model} object to estimate parameters
##'     in.
##' @slot multi_model FIXME.
##' @slot n_models Integer with the number of models.
##' @export
setClass(
    "SimInf_multi_model",
    slots = c(model       = "SimInf_model",
              multi_model = "SimInf_model",
              n_models    = "integer")
)
