## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2017  Stefan Engblom
## Copyright (C) 2015 - 2017  Stefan Widgren
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

library("SimInf")

## For debugging
sessionInfo()

m <- mparse(c("S -> b*S*I/(S+I+R) -> I", "I -> g*I -> R"),
            c("S", "I", "R"), list(b = 0.16, g = 0.077))

path <- tempdir()
package_skeleton(m, name = "SIR", path = path)

stopifnot(file.exists(file.path(path, "SIR", "DESCRIPTION")))
stopifnot(file.exists(file.path(path, "SIR", "NAMESPACE")))
stopifnot(file.exists(file.path(path, "SIR", "man", "SIR-class.Rd")))
stopifnot(file.exists(file.path(path, "SIR", "man", "SIR.Rd")))
stopifnot(file.exists(file.path(path, "SIR", "man", "run-methods.Rd")))
stopifnot(file.exists(file.path(path, "SIR", "R", "model.R")))
stopifnot(file.exists(file.path(path, "SIR", "src", "model.c")))
