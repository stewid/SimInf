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

##' plot,-method
##'
##'
##' @name plot-methods
##' @aliases plot plot-methods plot,SISe-method
##' @docType methods
##' @importFrom graphics plot
##' @export
setMethod("plot",
          signature(x = "SISe"),
          function(x, ...)
      {
          savepar <- par(mfrow = c(3,1),
                         mar = c(2,4,1,1),
                         oma = c(2,1,0,0))
          on.exit(par(savepar))

          I <- colSums(infected(x))
          S <- colSums(susceptible(x))

          plot(I / (S + I), type = "l", ylab = "Prevalence")
          plot(I, type = "l", ylab = "Infected")
          plot(S, t = "l", ylab = "Susceptible")

          title(xlab = "Day", outer = TRUE, line = 0)
      }
)
