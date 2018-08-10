## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2018  Stefan Engblom
## Copyright (C) 2015 - 2018  Stefan Widgren
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

## Check missing and invalid model argument
res <- tools::assertError(package_skeleton())
stopifnot(length(grep("Missing 'model' argument",
                      res[[1]]$message)) > 0)
res <- tools::assertError(package_skeleton(5))
stopifnot(length(grep("'model' argument is not a 'SimInf_model' object",
                      res[[1]]$message)) > 0)

## Chack package_skeleton
m <- mparse(transitions = c("S -> b*S*I/(S+I+R) -> I", "I -> g*I -> R"),
            compartments = c("S", "I", "R"),
            gdata = c(b = 0.16, g = 0.077),
            u0 = data.frame(S = 99, I = 1, R = 0),
            tspan = 1:10)

path <- tempdir()
package_skeleton(m, name = "SIR", path = path)

stopifnot(file.exists(file.path(path, "SIR", "DESCRIPTION")))
stopifnot(file.exists(file.path(path, "SIR", "NAMESPACE")))
stopifnot(file.exists(file.path(path, "SIR", "man", "SIR-class.Rd")))
stopifnot(file.exists(file.path(path, "SIR", "man", "SIR.Rd")))
stopifnot(file.exists(file.path(path, "SIR", "man", "run-methods.Rd")))
stopifnot(file.exists(file.path(path, "SIR", "R", "model.R")))
stopifnot(file.exists(file.path(path, "SIR", "src", "model.c")))

## Check that it fails if path exists
res <- tools::assertError(package_skeleton(m, name = "SIR", path = path))
stopifnot(length(grep("already exists",
                      res[[1]]$message)) > 0)

## Cleanup
unlink(path, recursive=TRUE)
