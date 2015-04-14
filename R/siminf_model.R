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

##' Class \code{"siminf_model"}
##'
##' Class to handle the siminf data model
##' @section Slots:
##' \describe{
##'   \item{G}{
##'     Sparse matrix (\eqn{Nt \times Nt}) of object class
##'     \code{"\linkS4class{dgCMatrix}"}.  A non-zeros entry in element
##'     \code{i} of column \code{j} indicates that propensity \code{i}
##'     needs to be recalculated if the transition \code{j} occurs.
##'   }
##'   \item{N}{
##'     Sparse matrix (\eqn{Nc \times Nt}) of object class
##'     \code{"\linkS4class{dgCMatrix}"}. Each column corresponds
##'     to a transition, and execution of transition \code{j} amounts to
##'     adding the \code{j}th column to the state vector.
##'   }
##'   \item{U}{
##'     The result matrix ((Nn * Nc) X length(tspan)). U(:,j) contains
##'     the state of the system at tspan(j).
##'   }
##'   \item{Nn}{
##'     Number of nodes.
##'   }
##'   \item{data}{
##'     A data vector passed as an argument to the propensities.
##'   }
##'   \item{sd}{
##'     Integer vector of length Nn. Each node can be assigned to a sub-domain.
##'   }
##'   \item{tspan}{
##'     A vector of increasing time points where the state of each node is
##'     to be returned.
##'   }
##'   \item{u0}{
##'      Initial state vector u0. Integer (Nc X Nn). Gives the initial
##'      number of individuals in each compartment in every node.
##'   }
##'   \item{events}{
##'     External events \code{"\linkS4class{external_events}"}
##'   }
##' }
##' @name siminf_model-class
##' @include external_events.r
##' @docType class
##' @keywords classes
##' @keywords methods
##' @export
##' @import Matrix
setClass("siminf_model",
         slots = c(G      = "dgCMatrix",
                   N      = "dgCMatrix",
                   U      = "matrix",
                   Nn     = "integer",
                   data   = "matrix",
                   sd     = "integer",
                   tspan  = "numeric",
                   u0     = "matrix",
                   events = "external_events"),
         validity = function(object) {
             errors <- character()

             ## Check Nn
             if (!identical(length(object@Nn), 1L)) {
               errors <- c(errors, "Wrong size of Nn.")
             }

             ## Check tspan.
             if (!is.double(object@tspan)) {
                 errors <- c(errors, "Input time-span must be a double vector.")
             } else if (any(length(object@tspan) < 2,
                            any(diff(object@tspan) <= 0))) {
                 errors <- c(errors, "Input time-span must be an increasing vector.")
             }

             ## Check u0.
             if (!identical(storage.mode(object@u0), "integer")) {
                 errors <- c(errors, "Initial state must be an integer matrix.")
             } else if (any(object@u0 < 0L)) {
                 errors <- c(errors, "Initial state has negative elements.")
             }

             ## Check U
             if (!identical(storage.mode(object@U), "integer")) {
                 errors <- c(errors, "Output state must be an integer matrix.")
             } else if (any(object@U < 0L)) {
                 errors <- c(errors, "Output state has negative elements.")
             }

             ## Check N.
             if (!all(is_wholenumber(object@N@x))) {
               stop("Stochiometric matrix must be an integer matrix.")
             }

             ## Check G.
             Nt <- dim(object@N)[2]
             if (!identical(dim(object@G), c(Nt, Nt))) {
                 errors <- c(errors, "Wrong size of dependency graph.")
             }

             ## Check sd.
             if (!identical(length(object@sd), object@Nn[1])) {
                 errors <- c(errors, "Wrong size of subdomain vector.")
             }

             ## Check data.
             if (!is.double(object@data)) {
                 errors <- c(errors, "Data matrix must be a double matrix.")
             }
             if (!identical(dim(object@data)[2], object@Nn[1])) {
                 errors <- c(errors, "Wrong size of data matrix.")
             }

             if (length(errors) == 0) TRUE else errors
         }
)

##' Create a siminf model
##'
##' @param G Sparse matrix (\eqn{Nt \times Nt}) of object class
##' \code{"\linkS4class{dgCMatrix}"}.  A non-zeros entry in element
##' \code{i} of column \code{j} indicates that propensity \code{i}
##' needs to be recalculated if the transition \code{j} occurs.
##' @param N Sparse matrix (\eqn{Nc \times Nt}) of object class
##' \code{"\linkS4class{dgCMatrix}"}. Each column corresponds to a
##' transition, and execution of transition \code{j} amounts to adding
##' the \code{j}th column to the state vector.
##' @param U The result matrix ((Nn * Nc) X length(tspan)). U(:,j)
##' contains the state of the system at tspan(j).
##' @param Nn Number of nodes.
##' @param data A data vector passed as an argument to the
##' propensities.
##' @param sd Integer vector of length Nn. Each node can be assigned
##' to a sub-domain.
##' @param tspan A vector of increasing time points where the state of
##' each node is to be returned.
##' @param u0 Initial state vector u0. Integer (Nc X Nn). Gives the
##' initial number of individuals in each compartment in every node.
##' @param events A \code{data.frame} with the scheduled events.
##' @param init A \code{data.frame} with the initial number of
##' individuals in each compartment in every node.
##' @param E Sparse matrix to handle external events, see
##' \code{\linkS4class{external_events}}.
##' @param S Sparse matrix to handle external events, see
##' \code{\linkS4class{external_events}}.
##' @return \linkS4class{siminf_model}
##' @export
siminf_model <- function(G,
                         N,
                         tspan,
                         events = NULL,
                         sd     = NULL,
                         data   = NULL,
                         U      = NULL,
                         Nn     = NULL,
                         u0     = NULL,
                         init   = NULL,
                         E      = NULL,
                         S      = NULL)
{
    ## Check initial state
    if (all(is.null(u0), is.null(init)))
        stop("Both u0 and init are NULL")
    if (all(!is.null(u0), !is.null(init)))
        stop("Both u0 and init are non NULL")

    ## Check u0
    if (!is.null(u0)) {
        if (!all(is.matrix(u0), is.numeric(u0)))
            stop("u0 must be an integer matrix")
        if (!is.integer(u0)) {
            if (!all(is_wholenumber(u0)))
                stop("u0 must be an integer matrix")
            storage.mode(u0) <- "integer"
        }
        if (is.null(Nn))
            Nn <- ncol(u0)
    }

    ## Check init
    if (!is.null(init)) {
        if (!is.data.frame(init))
            stop("init must be a data.frame")
        if (!("id" %in% names(init)))
            stop("init must contain the column id")
        if (!is.integer(init$id)) {
            if (!all(is_wholenumber(init$id)))
                stop("init$id must be an integer")
            init$id <- as.integer(init$id)
        }
        init <- init[order(init$id),]
        if (!identical(min(init$id), 0L))
            stop("init$id must be zero based")
        if (!identical(init$id, seq_len(max(init$id+1L))-1L))
            stop("init$id must be a sequence from 0 to max(init$id)-1")

        init$id <- NULL
        n.col <- ncol(init)
        if (is.null(Nn))
            Nn <- nrow(init)
        init <- t(data.matrix(init))
        attributes(init) <- NULL
        dim(init) <- c(n.col, Nn)
        u0 <- init
        storage.mode(u0) <- "integer"
    }

    ## Check Nn
    if (!is.integer(Nn)) {
        if (!all(is_wholenumber(Nn)))
            stop("Nn must be an integer")
        Nn <- as.integer(Nn)
    }
    if (!identical(Nn, ncol(u0)))
        stop("Nn must be equal to number of nodes")

    ## Check G
    if (class(G) == "dsCMatrix")
        G <- as(G, "dgCMatrix")

    ## Check data
    if (is.null(data))
        data <- matrix(rep(0, Nn), nrow = 1)

    ## Check U
    if (is.null(U)) {
        U <- matrix(nrow = 0, ncol = 0)
        storage.mode(U) <- "integer"
    } else {
        if (!is.integer(U)) {
            if (!all(is_wholenumber(U)))
                stop("U must be an integer")
            U <- as.integer(U)
        }

        if (!is.matrix(U)) {
            if (!identical(length(U), 0L))
                stop("U must be equal to 0 x 0 matrix")
            dim(U) <- c(0, 0)
        }
    }

    ## Check sd
    if (is.null(sd))
        sd <- rep(0L, Nn)

    ## Check events
    if (any(is.null(events), is.data.frame(events)))
        events <- external_events(E = E, S = S, events = events)

    return(new("siminf_model",
               G            = G,
               N            = N,
               U            = U,
               Nn           = Nn,
               data         = data,
               sd           = sd,
               tspan        = as.numeric(tspan),
               u0           = u0,
               events       = events))
}

##' Brief summary of \code{siminf_model}
##'
##' @aliases show,siminf_model-methods
##' @docType methods
##' @param object The siminf_model \code{object}
##' @return None (invisible 'NULL').
##' @keywords methods
##' @export
setMethod("show",
          signature(object = "siminf_model"),
          function (object)
          {
              cat("Epidemiological model:\n")
              cat(sprintf("G: %i x %i\n", dim(object@G)[1], dim(object@G)[2]))
              cat(sprintf("N: %i x %i\n", dim(object@N)[1], dim(object@N)[2]))
              cat(sprintf("U: %i x %i\n", dim(object@U)[1], dim(object@U)[2]))
              cat(sprintf("Nn: %i\n", object@Nn))
              cat(sprintf("data: %i x %i\n", dim(object@data)[1], dim(object@data)[2]))
              cat(sprintf("sd: %i x %i\n", dim(object@sd)[1], dim(object@sd)[2]))
              cat(sprintf("tspan: 1 x %i\n", length(object@tspan)))
              cat(sprintf("u0: %i x %i\n", dim(object@u0)[1], dim(object@u0)[2]))

              show(object@events)
          }
)
