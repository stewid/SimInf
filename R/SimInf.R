## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015  Stefan Engblom
## Copyright (C) 2015  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' A Framework for Stochastic Disease Spread Simulations
##'
##' @docType package
##' @name SimInf
##' @import methods
##' @useDynLib SimInf, .registration=TRUE
NULL

##' Unload hook function
##'
##' @param libpath A character string giving the complete path to the
##' package.
##' @keywords internal
.onUnload <- function (libpath)
{
    library.dynam.unload("SimInf", libpath)
}
