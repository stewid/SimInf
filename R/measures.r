## siminf, a framework for stochastic disease spread simulations
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

##' Susceptible
##'
##' Extracts the number of susceptible
##' @rdname susceptible-methods
##' @docType methods
##' @param model The \code{model} to extract the susceptible from
##' @param ... Additional arguments affecting the measure
##' @param age The age category to extract
##' ##' @keywords methods
##' @export
setGeneric("susceptible",
           function(model, ...) standardGeneric("susceptible"))

##' @rdname susceptible-methods
##' @include SISe.r
##' @export
setMethod("susceptible",
          signature("SISe"),
          function(model, ...) {
              as.matrix(model@U[seq(from = 1, to = dim(model@U)[1], by = 2), ])
          })

##' @rdname susceptible-methods
##' @include SISe3.r
##' @export
setMethod("susceptible",
          signature("SISe3"),
          function(model, age = c("age_1", "age_2", "age_3"), ...) {
              age <- match.arg(age)
              from <- switch(age,
                             age_1 = 1,
                             age_2 = 3,
                             age_3 = 5)
              to = dim(model@U)[1]
              as.matrix(model@U[seq(from = from, to = to, by = 6), ])
          })

##' Infected
##'
##' Extracts the number of infected
##' @rdname infected-methods
##' @docType methods
##' @param model The \code{model} to extract the infected from
##' @param ... Additional arguments affecting the measure
##' @param age The age category to extract
##' @keywords methods
##' @export
setGeneric("infected",
           function(model, ...) standardGeneric("infected"))

##' @rdname infected-methods
##' @include SISe.r
##' @export
setMethod("infected",
          signature("SISe"),
          function(model, ...) {
              as.matrix(model@U[seq(from = 2, to = dim(model@U)[1], by = 2), ])
          })

##' @rdname infected-methods
##' @include SISe3.r
##' @export
setMethod("infected",
          signature("SISe3"),
          function(model, age = c("age_1", "age_2", "age_3"), ...) {
              age <- match.arg(age)
              from <- switch(age,
                             age_1 = 2,
                             age_2 = 4,
                             age_3 = 6)
              to = dim(model@U)[1]
              as.matrix(model@U[seq(from = from, to = to, by = 6), ])
          })
