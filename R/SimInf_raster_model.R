## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
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

##' Create a \code{SimInf_raster_model}
##'
##' Create a \code{SimInf_raster_model} object.  It is a model where
##' the nodes are not fixed at one position but can move between cells
##' on a raster.
##' @template G-param
##' @template S-param
##' @param cell_S FIXME.
##' @template tspan-param
##' @template ldata-param
##' @template gdata-param
##' @template u0-param
##' @template v0-param
##' @template C_code-param
##' @export
SimInf_raster_model <- function(G,
                                S,
                                cell_S,
                                tspan,
                                ldata  = NULL,
                                gdata  = NULL,
                                u0     = NULL,
                                v0     = NULL,
                                C_code = NULL) {

    model <- SimInf_model(G      = G,
                          S      = S,
                          tspan  = tspan,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0,
                          C_code = C_code)

    methods::as(model, "SimInf_raster_model")
}
