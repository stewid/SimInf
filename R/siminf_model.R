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
##'   \item{S}{
##'     Sparse matrix (\eqn{Nc \times Nt}) of object class
##'     \code{"\linkS4class{dgCMatrix}"}. Each column corresponds
##'     to a transition, and execution of transition \code{j} amounts to
##'     adding the \code{j}th column to the state vector.
##'   }
##'   \item{U}{
##'     The compartment result matrix ((Nn * Nc) X length(tspan)).
##'     U(:,j) contains the state of the compartments at tspan(j).
##'   }
##'   \item{ldata}{
##'     A matrix with local data for nodes. The column ldata[, j] contains
##'     the local data vector for node #j. The local data vector is passed
##'     as an argument to the propensities and the post time step function.
##'   }
##'   \item{gdata}{
##'     A vector with global data for all nodes. The global data vector
##'     is passed as an argument to the propensities and the post time
##'     step function.
##'   }
##'   \item{sd}{
##'     Integer vector of length Nn. Each node can be assigned to a sub-domain.
##'   }
##'   \item{tspan}{
##'     A vector of increasing time points where the state of each node is
##'     to be returned.
##'   }
##'   \item{u0}{
##'     Initial compartment state vector u0. Integer (Nc X Nn). Gives the
##'     initial number of individuals in each compartment in every node.
##'   }
##'   \item{V}{
##'     The model parameter result matrix for parameters that are modelled
##'     ((Nn * Number of paramters) X length(tspan)).
##'     V(:,j) contains the state of the model parameters at tspan(j).
##'   }
##'   \item{v0}{
##'      Initial model state vector v0. Double (Number of parameters X Nn).
##'      Gives the initial value for model paramenters that are updated
##'      during the simulation in every node.
##'   }
##'   \item{events}{
##'     Scheduled events \code{"\linkS4class{scheduled_events}"}
##'   }
##' }
##' @name siminf_model-class
##' @include scheduled_events.R
##' @keywords classes
##' @keywords methods
##' @export
##' @import Matrix
setClass("siminf_model",
         slots = c(G      = "dgCMatrix",
                   S      = "dgCMatrix",
                   U      = "matrix",
                   ldata  = "matrix",
                   gdata  = "numeric",
                   sd     = "integer",
                   tspan  = "numeric",
                   u0     = "matrix",
                   V      = "matrix",
                   v0     = "matrix",
                   events = "scheduled_events"),
         validity = function(object) {
             errors <- character()

             ## Check tspan.
             if (!is.double(object@tspan)) {
                 errors <- c(errors, "Input time-span must be a double vector.")
             } else if (any(length(object@tspan) < 2,
                            any(diff(object@tspan) <= 0))) {
                 errors <- c(errors, "Input time-span must be an increasing vector.")
             }

             ## Check u0.
             if (!identical(storage.mode(object@u0), "integer")) {
                 errors <- c(errors, "Initial state 'u0' must be an integer matrix.")
             } else if (any(object@u0 < 0L)) {
                 errors <- c(errors, "Initial state 'u0' has negative elements.")
             }

             ## Check U.
             if (!identical(storage.mode(object@U), "integer")) {
                 errors <- c(errors, "Output state 'U' must be an integer matrix.")
             } else if (any(object@U < 0L)) {
                 errors <- c(errors, "Output state 'U' has negative elements.")
             }

             ## Check v0.
             if (!identical(storage.mode(object@v0), "double")) {
                 errors <- c(errors, "Initial model state 'v0' must be a double matrix.")
             }

             ## Check V.
             if (!identical(storage.mode(object@V), "double")) {
                 errors <- c(errors, "Output model state 'V' must be a double matrix.")
             }

             ## Check S.
             if (!all(is_wholenumber(object@S@x))) {
                 stop("Stochiometric matrix must be an integer matrix.")
             }

             ## Check G.
             Nt <- dim(object@S)[2]
             if (!identical(dim(object@G), c(Nt, Nt))) {
                 errors <- c(errors, "Wrong size of dependency graph.")
             }

             ## Check sd.
             Nn <- dim(object@u0)[2]
             if (!identical(length(object@sd), Nn)) {
                 errors <- c(errors, "Wrong size of subdomain vector.")
             }

             ## Check ldata.
             if (!is.double(object@ldata)) {
                 errors <- c(errors, "'ldata' matrix must be a double matrix.")
             }
             if (!identical(dim(object@ldata)[2], Nn)) {
                 errors <- c(errors, "Wrong size of 'ldata' matrix.")
             }

             ## Check gdata.
             if (!is.double(object@gdata)) {
                 errors <- c(errors, "'gdata' must be a double vector.")
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
##' @param S Sparse matrix (\eqn{Nc \times Nt}) of object class
##' \code{"\linkS4class{dgCMatrix}"}. Each column corresponds to a
##' transition, and execution of transition \code{j} amounts to adding
##' the \code{j}th column to the state vector.
##' @param U The result matrix ((Nn * Nc) X length(tspan)). U(:,j)
##' contains the state of the system at tspan(j).
##' @param ldata A matrix with local data for nodes. The column
##' ldata[, j] contains the local data vector for node #j. The local
##' data vector is passed as an argument to the propensities and the
##' post time step function.
##' @param gdata A vector with global data for all nodes. The global
##' data vector is passed as an argument to the propensities and the
##' post time step function.
##' @param sd Integer vector of length Nn. Each node can be assigned
##' to a sub-domain.
##' @param tspan A vector of increasing time points where the state of
##' each node is to be returned.
##' @param u0 Initial state vector u0. Integer (Nc X Nn). Gives the
##' initial number of individuals in each compartment in every node.
##' @param events A \code{data.frame} with the scheduled events.
##' @param init A \code{data.frame} with the initial number of
##' individuals in each compartment in every node.
##' @param V The model parameter result matrix for the parameters that
##' are modelled ((Nn * Number of paramters) X length(tspan)).  V(:,j)
##' contains the state of the model parameters at tspan(j).
##' @param v0 Initial model state vector v0. Double (Number of
##' parameters X Nn).  Gives the initial value for model paramenters
##' that are updated during the simulation in every node.
##' @param E Sparse matrix to handle scheduled events, see
##' \code{\linkS4class{scheduled_events}}.
##' @param N Sparse matrix to handle scheduled events, see
##' \code{\linkS4class{scheduled_events}}.
##' @return \linkS4class{siminf_model}
##' @export
siminf_model <- function(G,
                         S,
                         tspan,
                         events = NULL,
                         sd     = NULL,
                         ldata  = NULL,
                         gdata  = NULL,
                         U      = NULL,
                         u0     = NULL,
                         init   = NULL,
                         v0     = NULL,
                         V      = NULL,
                         E      = NULL,
                         N      = NULL)
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
        u0 <- t(data.matrix(init))
        attributes(u0) <- NULL
        dim(u0) <- c(ncol(init), nrow(init))
        storage.mode(u0) <- "integer"
    }

    ## Check G
    if (class(G) == "dsCMatrix")
        G <- as(G, "dgCMatrix")

    ## Check ldata
    if (is.null(ldata))
        ldata <- matrix(rep(0, ncol(u0)), nrow = 1)

    ## Check gdata
    if (is.null(gdata))
        gdata <- numeric(0)

    ## Check U
    if (is.null(U)) {
        U <- matrix(nrow = 0, ncol = 0)
        storage.mode(U) <- "integer"
    } else {
        if (!is.integer(U)) {
            if (!all(is_wholenumber(U)))
                stop("U must be an integer")
            storage.mode(U) <- "integer"
        }

        if (!is.matrix(U)) {
            if (!identical(length(U), 0L))
                stop("U must be equal to 0 x 0 matrix")
            dim(U) <- c(0, 0)
        }
    }

    ## Check v0
    if (is.null(v0)) {
        v0 <- matrix(nrow = 0, ncol = 0)
        storage.mode(v0) <- "double"
    } else {
        if (!all(is.matrix(v0), is.numeric(v0)))
            stop("v0 must be a numeric matrix")

        if (!identical(storage.mode(v0), "double"))
            storage.mode(v0) <- "double"
    }

    ## Check V
    if (is.null(V)) {
        V <- matrix(nrow = 0, ncol = 0)
        storage.mode(V) <- "double"
    } else {
        if (!is.numeric(V))
            stop("V must be numeric")

        if (!identical(storage.mode(V), "double"))
            storage.mode(V) <- "double"

        if (!is.matrix(V)) {
            if (!identical(length(V), 0L))
                stop("V must be equal to 0 x 0 matrix")
            dim(V) <- c(0, 0)
        }
    }

    ## Check sd
    if (is.null(sd))
        sd <- rep(0L, ncol(u0))

    ## Check events
    if (any(is.null(events), is.data.frame(events)))
        events <- scheduled_events(E = E, N = N, events = events)

    return(new("siminf_model",
               G      = G,
               S      = S,
               U      = U,
               ldata  = ldata,
               gdata  = gdata,
               sd     = sd,
               tspan  = as.numeric(tspan),
               u0     = u0,
               v0     = v0,
               V      = V,
               events = events))
}

##' Plot \code{\linkS4class{siminf_model}}
##'
##' @param x The \code{model} to plot
##' @param legend The character vector to appear in the legend.
##' @param t0 The first date of \code{x@@tspan} as a character string
##' in format 'yyyy-mm-dd'. Default is NULL which prints the x-axis
##' labels as the sequence 1:length(x@@tspan). If non-null, the labels
##' are converted to dates.
##' @param col The plotting color for each compartment. Default is
##' black.
##' @param lty The line type for each compartment. Default is the
##' sequence: 1=solid, 2=dashed, 3=dotted, 4=dotdash, 5=longdash,
##' 6=twodash.
##' @param ... Additional arguments affecting the plot produced.
##' @name plot-methods
##' @aliases plot plot-methods plot,siminf_model-method
##' @docType methods
##' @importFrom graphics axis
##' @importFrom graphics legend
##' @importFrom graphics lines
##' @importFrom graphics par
##' @importFrom graphics plot
##' @importFrom graphics title
##' @export
##' @examples
##' ## Create a 'SISe' demo model with 1 node and initialize
##' ## it to run over 1000 days.
##' model <- demo_model(nodes = 1, days = 1000, model = "SISe")
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Plot the proportion susceptible and infected individuals
##' plot(result)
setMethod("plot",
          signature(x = "siminf_model"),
          function(x, legend, t0 = NULL, col = NULL, lty = NULL, ...)
      {
          if (identical(dim(x@U), c(0L, 0L)))
              stop("Please run the model first, the 'U' matrix is empty")

          savepar <- par(mar = c(2,4,1,1), oma = c(4,1,0,0), xpd = TRUE)
          on.exit(par(savepar))

          ## Create matrix where each row is the sum of individuals in
          ## that state
          m <- do.call(rbind, lapply(seq_len(dim(x@S)[1]), function(from) {
              i <- seq(from = from, to = dim(x@U)[1], by = dim(x@S)[1])
              colSums(as.matrix(x@U[i, , drop = FALSE]))
          }))

          ## Calculate proportion
          m <- m / colSums(m)

          ## Default color is black
          if (is.null(col))
              col <- rep("black", dim(x@S)[1])

          ## Default line type
          if (is.null(lty))
              lty <- seq_len(dim(m)[1])

          ## Plot
          if (is.null(t0)) {
              plot(m[1,], type = "l", ylab = "Proportion", ylim = c(0, max(m)),
                   col = col[1], lty = lty[1], ...)
          } else {
              plot(m[1,], type = "l", ylab = "Proportion", ylim = c(0, max(m)),
                   col = col[1], lty = lty[1], xaxt = "n")
          }
          title(xlab = "Day", outer = TRUE, line = 0)
          if (!is.null(t0)) {
              axis(side = 1, at = seq_len(dim(m)[2]),
                   labels = as.Date(x@tspan - min(x@tspan),
                       origin = as.Date(t0)))
          }
          for (i in seq_len(dim(m)[1])[-1]) {
              lines(m[i, ], type = "l", lty = lty[i], col = col[i], ...)
          }

          ## Add legend below plot
          par(fig = c(0, 1, 0, 1),
              oma = c(0, 0, 0, 0),
              mar = c(0, 0, 0, 0), new = TRUE)
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          graphics::legend("bottom", inset = c(0, 0),
                           lty = lty, col = col, bty = "n",
                           horiz = TRUE, legend = legend)
      }
)

##' Brief summary of \code{siminf_model}
##'
##' @aliases show,siminf_model-methods
##' @docType methods
##' @param object The siminf_model \code{object}
##' @return None (invisible 'NULL').
##' @keywords methods
##' @export
##' @examples
##' ## Create a 'SISe' demo model with 1 node and initialize
##' ## it to run over 1000 days.
##' model <- demo_model(nodes = 1, days = 1000, model = "SISe")
##'
##' ## Brief summary of the model
##' model
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Brief summary of the result.
##' result
setMethod("show",
          signature(object = "siminf_model"),
          function (object)
          {
              cat("Epidemiological model:\n")
              cat(sprintf("G: %i x %i\n", dim(object@G)[1], dim(object@G)[2]))
              cat(sprintf("S: %i x %i\n", dim(object@S)[1], dim(object@S)[2]))
              cat(sprintf("U: %i x %i\n", dim(object@U)[1], dim(object@U)[2]))
              cat(sprintf("V: %i x %i\n", dim(object@V)[1], dim(object@V)[2]))
              cat(sprintf("ldata: %i x %i\n", dim(object@ldata)[1], dim(object@ldata)[2]))
              cat(sprintf("gdata: 1 x %i\n", length(object@gdata)))
              cat(sprintf("sd: %i x %i\n", dim(object@sd)[1], dim(object@sd)[2]))
              cat(sprintf("tspan: 1 x %i\n", length(object@tspan)))
              cat(sprintf("u0: %i x %i\n", dim(object@u0)[1], dim(object@u0)[2]))
              cat(sprintf("v0: %i x %i\n", dim(object@v0)[1], dim(object@v0)[2]))

              show(object@events)
          }
)
