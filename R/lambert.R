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

##' Lambert W0 function
##'
##' W0(x) is the principal branch of the solution of the function
##' defined by \eqn{We^W = x}{W * exp(W)} for \eqn{x >= -1/e}. The
##' value is calculated using GNU Scientific Library (GSL).
##' @param x numeric vector of values.
##' @references
##'
##' GNU Scientific Library <https://www.gnu.org/software/gsl/>
##' @export
##' @examples
##' ## Should equal 1, as 1 * exp(1) = e.
##' lambertW0(exp(1))
##'
##' ## Should equal 0, as 0 * exp(0) = 0.
##' lambertW0(0)
##'
##' ## Should equal -1.
##' lambertW0(-exp(-1))
lambertW0 <- function(x) {
    .Call(SimInf_lambertW0, x)
}
