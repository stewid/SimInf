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

##' Class \code{"SimInf_mparse"}
##'
##' Class to handle the SimInf mparse data
##' @section Slots:
##' \describe{
##'   \item{C_code}{
##'     Character vector with model C code.
##'   }
##'   \item{G}{
##'     Dependency graph that indicates the transition rates that need
##'     to be updated after a given state transition has occured.
##'     A non-zero entry in element \code{G[i, i]} indicates that transition
##'     rate \code{i} needs to be recalculated if the state transition
##'     \code{j} occurs. Sparse matrix (\eqn{Nt \times Nt}) of object class
##'     \code{"\linkS4class{dgCMatrix}"}.
##'   }
##'   \item{S}{
##'     Each column corresponds to a state transition, and execution
##'     of state transition \code{j} amounts to adding the \code{S[,
##'     j]} column to the state vector \code{u[, i]} of node \emph{i}
##'     where the transition occurred. Sparse matrix (\eqn{Nc \times
##'     Nt}) of object class \code{"\linkS4class{dgCMatrix}"}.
##'   }
##' }
##' @keywords methods
##' @export
##' @import Matrix
setClass("SimInf_mparse",
         slots = c(C_code = "character",
                   G      = "dgCMatrix",
                   S      = "dgCMatrix"),
         validity = function(object) {
             errors <- character()

             ## Check C code
             if (nchar(paste0(object@C_code, collapse = "\n")) == 0L) {
                 errors <- c(errors, "'C_code' is empty.")
             }

             ## Check S.
             if (identical(dim(object@S), c(0L, 0L))) {
                 errors <- c(errors, "'S' is empty.")
             } else if (!all(is_wholenumber(object@S@x))) {
                 errors <- c(errors,
                             "'S' matrix must be an integer matrix.")
             }

             ## Check G.
             Nt <- dim(object@S)[2]
             if (!identical(dim(object@G), c(Nt, Nt))) {
                 errors <- c(errors,
                             "Wrong size of dependency graph 'G'.")
             }

             if (length(errors) == 0) TRUE else errors
         }
)
