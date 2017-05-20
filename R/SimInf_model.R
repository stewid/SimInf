## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
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

##' Class \code{"SimInf_model"}
##'
##' Class to handle the siminf data model
##' @section Slots:
##' \describe{
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
##'   \item{U}{
##'     The result matrix with the number of individuals in each
##'     compartment in every node. \code{U[, j]} contains the number
##'     of individuals in each compartment at
##'     \code{tspan[j]}. \code{U[1:Nc, j]} contains the number of
##'     individuals in node 1 at \code{tspan[j]}. \code{U[(Nc + 1):(2
##'     * Nc), j]} contains the number of individuals in node 2 at
##'     \code{tspan[j]} etc. Integer matrix (\eqn{N_n N_c \times}
##'     \code{length(tspan)}).
##'   }
##'   \item{U_sparse}{
##'     If the model was run to write the solution to a sparse matrix
##'     (\code{dgCMatrix}) the \code{U_sparse} contains the data and
##'     \code{U} is empty. The layout of the data in \code{U_sparse}
##'     is identical to \code{U}. Please note that \code{U_sparse}
##'     is numeric and \code{U} is integer.
##'   }
##'   \item{V}{
##'     The result matrix for the real-valued continuous
##'     state. \code{V[, j]} contains the real-valued state of the
##'     system at \code{tspan[j]}. Numeric matrix
##'     (\eqn{N_n}\code{dim(ldata)[1]} \eqn{\times}
##'     \code{length(tspan)}).
##'   }
##'   \item{V_sparse}{
##'     If the model was run to write the solution to a sparse matrix
##'     (\code{dgCMatrix}) the \code{V_sparse} contains the data and
##'     \code{V} is empty. The layout of the data in \code{V_sparse}
##'     is identical to \code{V}.
##'   }
##'   \item{ldata}{
##'     A matrix with local data for the nodes. The column \code{ldata[, j]}
##'     contains the local data vector for the node \code{j}. The local
##'     data vector is passed as an argument to the transition rate
##'     functions and the post time step function.
##'   }
##'   \item{gdata}{
##'     A numeric vector with global data that is common to all nodes.
##'     The global data vector is passed as an argument to the
##'     transition rate functions and the post time step function.
##'   }
##'   \item{tspan}{
##'     A vector of increasing time points where the state of each node is
##'     to be returned.
##'   }
##'   \item{u0}{
##'     The initial state vector (\eqn{N_c \times N_n}) with
##'     the number of individuals in each compartment in every node.
##'   }
##'   \item{v0}{
##'      The initial value for the real-valued continuous state.
##'      Numeric matrix (\code{dim(ldata)[1]} \eqn{\times N_n}).
##'   }
##'   \item{events}{
##'     Scheduled events \code{"\linkS4class{SimInf_events}"}
##'   }
##'   \item{C_code}{
##'     Character vector with optional model C code. If non-empty, the
##'     C code is written to a temporary C-file when the \code{run}
##'     method is called.  The temporary C-file is compiled and the
##'     resulting DLL is dynamically loaded. The DLL is unloaded and
##'     the temporary files are removed after running the model.
##'   }
##' }
##' @include SimInf_events.R
##' @keywords methods
##' @export
##' @import Matrix
setClass("SimInf_model",
         slots = c(G        = "dgCMatrix",
                   S        = "dgCMatrix",
                   U        = "matrix",
                   U_sparse = "dgCMatrix",
                   ldata    = "matrix",
                   gdata    = "numeric",
                   tspan    = "numeric",
                   u0       = "matrix",
                   V        = "matrix",
                   V_sparse = "dgCMatrix",
                   v0       = "matrix",
                   events   = "SimInf_events",
                   C_code   = "character"),
         validity = function(object) {
             ## Check events
             errors <- validObject(object@events)
             if (identical(errors, TRUE))
                 errors <- character()

             ## Check tspan.
             if (!is.double(object@tspan)) {
                 errors <- c(errors, "Input time-span must be a double vector.")
             } else if (any(length(object@tspan) < 2,
                            any(diff(object@tspan) <= 0),
                            any(is.na(object@tspan)))) {
                 errors <- c(errors,
                             "Input time-span must be an increasing vector.")
             }

             ## Check u0.
             if (!identical(storage.mode(object@u0), "integer")) {
                 errors <- c(errors,
                             "Initial state 'u0' must be an integer matrix.")
             } else if (any(object@u0 < 0L)) {
                 errors <- c(errors,
                             "Initial state 'u0' has negative elements.")
             }

             ## Check U.
             if (!identical(storage.mode(object@U), "integer")) {
                 errors <- c(errors,
                             "Output state 'U' must be an integer matrix.")
             } else if (any(object@U < 0L)) {
                 errors <- c(errors,
                             "Output state 'U' has negative elements.")
             }

             ## Check v0.
             if (!identical(storage.mode(object@v0), "double")) {
                 errors <- c(errors,
                             "Initial model state 'v0' must be a double matrix.")
             }

             ## Check V.
             if (!identical(storage.mode(object@V), "double")) {
                 errors <- c(errors,
                             "Output model state 'V' must be a double matrix.")
             }

             ## Check S.
             if (!all(is_wholenumber(object@S@x))) {
                 errors <- c(errors,
                             "'S' matrix must be an integer matrix.")
             }

             ## Check G.
             Nt <- dim(object@S)[2]
             if (!identical(dim(object@G), c(Nt, Nt))) {
                 errors <- c(errors,
                             "Wrong size of dependency graph.")
             }

             ## Check ldata.
             if (!is.double(object@ldata)) {
                 errors <- c(errors,
                             "'ldata' matrix must be a double matrix.")
             }
             Nn <- dim(object@u0)[2]
             if (!identical(dim(object@ldata)[2], Nn)) {
                 errors <- c(errors,
                             "Wrong size of 'ldata' matrix.")
             }

             ## Check gdata.
             if (!is.double(object@gdata)) {
                 errors <- c(errors,
                             "'gdata' must be a double vector.")
             }

             if (length(errors) == 0) TRUE else errors
         }
)

##' Create a \code{SimInf_model}
##'
##' @param G Dependency graph that indicates the transition rates that
##'     need to be updated after a given state transition has occured.
##'     A non-zero entry in element \code{G[i, i]} indicates that
##'     transition rate \code{i} needs to be recalculated if the state
##'     transition \code{j} occurs. Sparse matrix (\eqn{Nt \times Nt})
##'     of object class \code{"\linkS4class{dgCMatrix}"}.
##' @param S Each column corresponds to a transition, and execution of
##'     state transition \code{j} amounts to adding the \code{S[, j]}
##'     to the state vector of the node where the state transition
##'     occurred.  Sparse matrix (\eqn{Nc \times Nt}) of object class
##'     \code{"\linkS4class{dgCMatrix}"}.
##' @param U The result matrix with the number of individuals in each
##'     disease state in every node (\eqn{N_n N_c \times}
##'     \code{length(tspan)}).  \code{U[, j]} contains the number of
##'     individuals in each disease state at
##'     \code{tspan[j]}. \code{U[1:Nc, j]} contains the state of node
##'     \code{1} at \code{tspan[j]}. \code{U[(Nc + 1):(2 * Nc), j]}
##'     contains the state of node \code{2} at \code{tspan[j]} etc.
##' @param ldata A matrix with local data for the nodes. The column
##'     \code{ldata[, j]} contains the local data vector for the node
##'     \code{j}. The local data vector is passed as an argument to
##'     the transition rate functions and the post time step function.
##' @param gdata A numeric vector with global data that is common to
##'     all nodes. The global data vector is passed as an argument to
##'     the transition rate functions and the post time step function.
##' @template tspan-param
##' @param u0 The initial state vector. Either a matrix (\eqn{N_c
##'     \times N_n}) or a a \code{data.frame} with the number of
##'     individuals in each compartment in every node.
##' @param events A \code{data.frame} with the scheduled events.
##' @param V The result matrix for the real-valued continous
##'     compartment state (\eqn{N_n}\code{dim(ldata)[1]} \eqn{\times}
##'     \code{length(tspan)}).  \code{V[, j]} contains the real-valued
##'     state of the system at \code{tspan[j]}.
##' @param v0 The initial continuous state vector in every node.
##'     (\code{dim(ldata)[1]} \eqn{N_N \times}). The continuous state
##'     vector is updated by the specific model during the simulation
##'     in the post time step function.
##' @param E Sparse matrix to handle scheduled events, see
##'     \code{\linkS4class{SimInf_events}}.
##' @param N Sparse matrix to handle scheduled events, see
##'     \code{\linkS4class{SimInf_events}}.
##' @param C_code Character vector with optional model C code. If
##'     non-empty, the C code is written to a temporary C-file when
##'     the \code{run} method is called.  The temporary C-file is
##'     compiled and the resulting DLL is dynamically loaded. The DLL
##'     is unloaded and the temporary files are removed after running
##'     the model.
##' @return \linkS4class{SimInf_model}
##' @export
SimInf_model <- function(G,
                         S,
                         tspan,
                         events = NULL,
                         ldata  = NULL,
                         gdata  = NULL,
                         U      = NULL,
                         u0     = NULL,
                         v0     = NULL,
                         V      = NULL,
                         E      = NULL,
                         N      = NULL,
                         C_code = NULL)
{
    ## Check u0
    if (is.null(u0))
        stop("'u0' is NULL")
    if (is.data.frame(u0)) {
        n_col <- ncol(u0)
        n_row <- nrow(u0)
        lbl <- colnames(u0)
        u0 <- t(data.matrix(u0))
        attributes(u0) <- NULL
        dim(u0) <- c(n_col, n_row)
        rownames(u0) <- lbl
    }
    if (!all(is.matrix(u0), is.numeric(u0)))
        stop("u0 must be an integer matrix")
    if (!is.integer(u0)) {
        if (!all(is_wholenumber(u0)))
            stop("u0 must be an integer matrix")
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

    ## Check tspan
    if (is(tspan, "Date")) {
        ## Coerce the date vector to a numeric vector as days, where
        ## tspan[1] becomes the day of the year of the first year of
        ## the tspan date vector. The dates are added as names to the
        ## numeric vector.
        t0 <- as.numeric(as.Date(format(tspan[1], "%Y-01-01"))) - 1
        tspan_lbl <- format(tspan, "%Y-%m-%d")
        tspan <- as.numeric(tspan) - t0
        names(tspan) <- tspan_lbl
    } else {
        t0 <- NULL
    }
    storage.mode(tspan) <- "double"

    ## Check events
    if (!any(is.null(events), is.data.frame(events)))
        stop("'events' must be NULL or a data.frame")
    events <- SimInf_events(E = E, N = N, events = events, t0 = t0)

    ## Check C code
    if (is.null(C_code))
        C_code <- character(0)

    return(new("SimInf_model",
               G      = G,
               S      = S,
               U      = U,
               ldata  = ldata,
               gdata  = gdata,
               tspan  = tspan,
               u0     = u0,
               v0     = v0,
               V      = V,
               events = events,
               C_code = C_code))
}

## Internal function to calculate prevalence from U
calc_prevalence <- function(model = NULL, numerator = NULL,
                            denominator = NULL,
                            type = c("pop", "bnp", "wnp"), i = NULL)
{
    numerator <- extract_U(model, numerator, i)
    denominator <- extract_U(model, denominator, i)

    type <- match.arg(type)
    if (identical(type, "pop")) {
        numerator <- colSums(numerator)
        denominator <- colSums(denominator)
    } else if (identical(type, "bnp")) {
        numerator <- colSums(numerator > 0)
        ## Include only nodes with individuals
        denominator <- colSums(denominator > 0)
    }

    numerator / denominator
}

## Internal function to extract compartments from U
extract_U <- function(model = NULL, compartments = NULL, i = NULL) {
    if (is.null(model))
        stop("Missing 'model' argument")
    if (is.null(compartments))
        stop("Missing 'compartments' argument")
    if (!(all(compartments %in% rownames(model@S))))
        stop("'compartments' must exist in the model")
    compartments <- match(compartments, rownames(model@S))
    if (identical(dim(model@U), c(0L, 0L)))
        stop("Please run the model first, the 'U' matrix is empty")

    result <- NULL
    for (k in compartments) {
        ii <- seq(from = k, to = dim(model@U)[1], by = Nc(model))
        if (!is.null(i))
            ii <- ii[i]

        if (is.null(result)) {
            result <- as.matrix(model@U[ii, , drop = FALSE])
        } else {
            result <- result + as.matrix(model@U[ii, , drop = FALSE])
        }
    }

    rownames(result) <- NULL
    colnames(result) <- NULL
    result
}

##' @rdname U-methods
##' @export
setMethod("U",
          signature("SimInf_model"),
          function(model) {
              d <- dim(model@U)
              if (identical(d, c(0L, 0L))) {
                  d <- dim(model@U_sparse)
                  if (identical(d, c(0L, 0L)))
                      stop("Please run the model first, the 'U' matrix is empty")
                  return(model@U_sparse)
              }
              model@U
          }
)

##' @rdname U_set-methods
##' @export
setMethod("U<-",
          signature("SimInf_model"),
          function(model, value) {
              if (!is.null(value)) {
                  if (!is(value, "dgCMatrix"))
                      value <- as(value, "dgCMatrix")

                  d <- c(Nn(model) * Nc(model), length(model@tspan))
                  if (!identical(dim(value), d))
                      stop("Wrong dimension of 'value'")

                  ## Clear dense result matrix
                  u <- matrix(nrow = 0, ncol = 0)
                  storage.mode(u) <- "integer"
                  model@U = u

                  model@U_sparse = value
              } else {
                  ## Clear sparse result matrix
                  model@U_sparse <- as(sparseMatrix(numeric(0), numeric(0),
                                                    dims = c(0, 0)),
                                       "dgCMatrix")
              }
              model
          }
)

##' @rdname V-methods
##' @export
setMethod("V",
          signature("SimInf_model"),
          function(model) {
              d <- dim(model@V)
              if (identical(d, c(0L, 0L))) {
                  d <- dim(model@V_sparse)
                  if (identical(d, c(0L, 0L)))
                      stop("Please run the model first, the 'V' matrix is empty")
                  return(model@V_sparse)
              }
              model@V
          }
)

##' @rdname V_set-methods
##' @export
setMethod("V<-",
          signature("SimInf_model"),
          function(model, value) {
              if (!is.null(value)) {
                  if (!is(value, "dgCMatrix"))
                      value <- as(value, "dgCMatrix")

                  d <- c(Nn(model) * Nd(model), length(model@tspan))
                  if (!identical(dim(value), d))
                      stop("Wrong dimension of 'value'")

                  ## Clear dense result matrix
                  v <- matrix(nrow = 0, ncol = 0)
                  storage.mode(v) <- "double"
                  model@V <- v

                  model@V_sparse = value
              } else {
                  ## Clear sparse result matrix
                  model@V_sparse <- as(sparseMatrix(numeric(0), numeric(0),
                                                    dims = c(0, 0)),
                                       "dgCMatrix")
              }
              model
          }
)

## Number of nodes
Nn <- function(model) {
    dim(model@u0)[2]
}

## Number of compartments
Nc <- function(model) {
    dim(model@S)[1]
}

## Number of transitions
Nt <- function(model) {
    dim(model@G)[1]
}

## Number of continuous state variables
Nd <- function(model) {
    dim(model@v0)[1]
}

##' @rdname run-methods
##' @export
setMethod("run",
          signature(model = "SimInf_model"),
          function(model, threads, seed)
          {
              ## Check that SimInf_model contains all data structures
              ## required by the siminf solver and that they make sense
              validObject(model);

              if (nchar(paste0(model@C_code, collapse = "\n"))) {
                  ## Write the C code to a temporary file
                  filename <- tempfile("SimInf-")
                  on.exit(unlink(paste0(filename,
                                        c(".c", ".o", .Platform$dynlib.ex))))
                  writeLines(model@C_code, con = paste0(filename, ".c"))

                  ## Include directive for "SimInf.h"
                  include <- system.file("include", package = "SimInf")
                  Sys.setenv(PKG_CPPFLAGS=sprintf("-I%s", shQuote(include)))

                  ## Compile the model C code using the running version of R.
                  wd <- setwd(dirname(filename))
                  cmd <- paste(shQuote(file.path(R.home(component="bin"), "R")),
                               "CMD SHLIB",
                               shQuote(paste0(basename(filename), ".c")))
                  compiled <- system(cmd, intern = TRUE)
                  setwd(wd)

                  ## Load DLL
                  lib <- paste0(filename, .Platform$dynlib.ext)
                  if (!file.exists(lib))
                      stop(compiled)
                  dll <- dyn.load(lib)
                  on.exit(dyn.unload(lib), add = TRUE)

                  ## Create expression to parse
                  expr <- ".Call(dll$SimInf_model_run, model, threads, seed)"
              } else {
                  ## The model name
                  name <- as.character(class(model))

                  ## The model C run function
                  run_fn <- paste0(name, "_run")

                  ## Create expression to parse
                  expr <- ".Call(run_fn, model, threads, seed, PACKAGE = 'SimInf')"
              }

              ## Run model
              eval(parse(text = expr))
          }
)

## Create a matrix where each column contains the individuals in one
## state.
by_compartment <- function(model) {
    if (identical(dim(model@U), c(0L, 0L)))
        stop("Please run the model first, the 'U' matrix is empty")

    m <- do.call(cbind, lapply(seq_len(dim(model@S)[1]), function(from) {
        i <- seq(from = from, to = dim(model@U)[1], by = dim(model@S)[1])
        as.numeric(model@U[i, ])
    }))

    colnames(m) <- rownames(model@S)
    m
}

##' Box Plots
##'
##' Produce box-and-whisker plot(s) of the number of individuals in
##' each compartment.
##' @param x The \code{model} to plot
##' @param ... Additional arguments affecting the plot produced.
##' @name boxplot-methods
##' @aliases boxplot boxplot-methods boxplot,SimInf_model-method
##' @export
setMethod("boxplot",
          signature(x = "SimInf_model"),
          function(x, ...)
          {
              graphics::boxplot(by_compartment(x), ...)
          }
)

##' Scatterplot Matrices
##'
##' A matrix of scatterplots with the number of individuals is
##' produced. The \code{ij}th scatterplot contains \code{x[,i]}
##' plotted against \code{x[,j]}.
##' @param x The \code{model} to plot
##' @param ... Additional arguments affecting the plot produced.
##' @name pairs-methods
##' @aliases pairs pairs-methods pairs,SimInf_model-method
##' @export
setMethod("pairs",
          signature(x = "SimInf_model"),
          function(x, ...)
          {
              graphics::pairs(by_compartment(x), ...)
          }
)

##' Plot \code{\linkS4class{SimInf_model}}
##'
##' @param x The \code{model} to plot
##' @param legend The character vector to appear in the
##'     legend. Default is to use the names of the compartments.
##' @param col The plotting color for each compartment. Default is
##'     black.
##' @param lty The line type for each compartment. Default is the
##'     sequence: 1=solid, 2=dashed, 3=dotted, 4=dotdash, 5=longdash,
##'     6=twodash.
##' @param lwd The line width for each compartment. Default is 2.
##' @param N if \code{TRUE}, the average number of individuals in each
##'     compartment, else the proportion of individuals in each
##'     compartment.  Default is \code{FALSE}.
##' @param compartments Character vector with the compartments in the
##'     model to include in the plot. Default is \code{NULL}
##'     i.e. include all compartments in the model.
##' @param spaghetti Plot one line for each node. Default is
##'     \code{FALSE}.
##' @param ... Additional arguments affecting the plot produced.
##' @name plot-methods
##' @aliases plot plot-methods plot,SimInf_model-method
##' @importFrom graphics legend
##' @importFrom graphics lines
##' @importFrom graphics par
##' @importFrom graphics plot
##' @importFrom graphics title
##' @export
##' @examples
##' ## Create an 'SIR' model with 10 nodes and initialise
##' ## it to run over 100 days.
##' model <- SIR(u0 = data.frame(S = rep(99, 10),
##'                              I = rep(1, 10),
##'                              R = rep(0, 10)),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Run the model and save the result.
##' result <- run(model)
##'
##' ## Plot the proportion of susceptible, infected and recovered
##' ## individuals.
##' plot(result)
##'
##' ## Plot the number of susceptible, infected and recovered
##' ## individuals.
##' plot(result, N = TRUE)
##'
##' ## Plot the number of infected individuals.
##' plot(result, compartments = "I", N = TRUE)
setMethod("plot",
          signature(x = "SimInf_model"),
          function(x, legend = NULL, col = NULL, lty = NULL, lwd = 2,
                   N = FALSE, compartments = NULL, spaghetti = FALSE,
                   ...)
          {
              if (identical(dim(x@U), c(0L, 0L)))
                  stop("Please run the model first, the 'U' matrix is empty")

              ## Determine the compartments to include in the plot
              if (is.null(compartments)) {
                  compartments <- seq_len(Nc(x))
              } else {
                  if (!(all(compartments %in% rownames(x@S))))
                      stop("'compartments' must exist in the model")
                  compartments <- match(compartments, rownames(x@S))
              }

              savepar <- graphics::par(mar = c(2,4,1,1), oma = c(4,1,0,0),
                                       xpd = TRUE)
              on.exit(graphics::par(savepar))

              ## Create a matrix with one row for each line in the
              ## plot.
              if (identical(spaghetti, TRUE)) {
                  i <- sort(as.numeric(sapply(compartments, "+",
                  (seq_len(Nn(x)) - 1) * Nc(x))))
                  m <- x@U[i, , drop = FALSE]

                  if (!identical(N, TRUE)) {
                      ## Calculate the proportion of individuals in
                      ## each compartment within each node.
                      for (i in (seq_len(Nn(x)) - 1)) {
                          n <- colSums(x@U[i * Nc(x) + seq_len(Nc(x)), ,
                                           drop = FALSE])
                          j <- i * length(compartments) +
                              seq_len(length(compartments))
                          m[j, ] <- m[j, , drop = FALSE] / n
                      }
                  }
              } else {
                  m <- matrix(0, nrow = length(compartments),
                              ncol = length(x@tspan))
                  for (i in seq_len(length(compartments))) {
                      j <- seq(from = compartments[i], to = dim(x@U)[1],
                               by = Nc(x))
                      if (N) {
                          m[i, ] <- colMeans(as.matrix(x@U[j, , drop = FALSE]))
                      } else {
                          m[i, ] <- colSums(as.matrix(x@U[j, , drop = FALSE]))
                      }
                  }

                  ## Calculate proportion
                  if (!identical(N, TRUE))
                      m <- apply(m, 2, function(x) x / sum(x))
              }

              ## Default line type
              if (is.null(lty)) {
                  if (is.null(col)) {
                      lty <- seq_len(length(compartments))
                  } else {
                      lty <- rep(1, length(compartments))
                  }
              }
              lty <- rep(lty, length.out = dim(m)[1])

              ## Default color is black
              if (is.null(col)) {
                  col <- rep("black", length(compartments))
              } else {
                  col <- col[compartments]
              }
              col <- rep(col, length.out = dim(m)[1])

              ## Settings for the y-axis.
              if (identical(N, TRUE)) {
                  ylab <- "N"
                  ylim <- c(0, max(m))
              } else {
                  ylab <- "Proportion"
                  ylim <- c(0, 1)
              }

              ## Settings for the x-axis
              if (is.null(names(x@tspan))) {
                  xx <- x@tspan
                  xlab <- "Time"
              } else {
                  xx <- as.Date(names(x@tspan))
                  xlab <- "Date"
              }

              ## Plot first line to get a new plot window
              graphics::plot(x = xx, y = m[1, ], type = "l",
                             ylab = ylab, ylim = ylim, col = col[1],
                             lty = lty[1], lwd = lwd, ...)
              graphics::title(xlab = xlab, outer = TRUE, line = 0)

              ## Add the rest of the lines to the plot
              for (i in seq_len(dim(m)[1])[-1]) {
                  graphics::lines(x = xx, y = m[i, ], type = "l",
                                  lty = lty[i], col = col[i], lwd = lwd,
                                  ...)
              }

              ## Add the legend below plot. The default legend is the
              ## names of the compartments.
              if (is.null(legend))
                  legend <- rownames(x@S)[compartments]
              graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                            mar = c(0, 0, 0, 0), new = TRUE)
              graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
              graphics::legend("bottom", inset = c(0, 0),
                               lty = lty[seq_len(length(compartments))],
                               col = col[seq_len(length(compartments))],
                               bty = "n", horiz = TRUE,
                               legend = legend, lwd = lwd)
          }
)

##' @keywords internal
show_U <- function(object) {
    d <- dim(object@U)
    if (identical(d, c(0L, 0L))) {
        d <- dim(object@U_sparse)
        if (identical(d, c(0L, 0L))) {
            cat("U: 0 x 0\n")
        } else {
            cat(sprintf("U: %i x %i (sparse)\n", d[1], d[2]))
        }
    } else {
        cat(sprintf("U: %i x %i\n", d[1], d[2]))
    }
}

##' @keywords internal
show_V <- function(object) {
    d <- dim(object@V)
    if (identical(d, c(0L, 0L))) {
        d <- dim(object@V_sparse)
        if (identical(d, c(0L, 0L))) {
            cat("V: 0 x 0\n")
        } else {
            cat(sprintf("V: %i x %i (sparse)\n", d[1], d[2]))
        }
    } else {
        cat(sprintf("V: %i x %i\n", d[1], d[2]))
    }
}

##' Brief summary of \code{SimInf_model}
##'
##' @aliases show,SimInf_model-methods
##' @param object The SimInf_model \code{object}
##' @return None (invisible 'NULL').
##' @keywords methods
##' @export
##' @examples
##' ## Create an 'SIR' model with 10 nodes and initialise
##' ## it to run over 100 days.
##' model <- SIR(u0 = data.frame(S = rep(99, 10),
##'                              I = rep(1, 10),
##'                              R = rep(0, 10)),
##'              tspan = 1:100,
##'              beta = 0.16,
##'              gamma = 0.077)
##'
##' ## Brief summary of the model
##' model
##'
##' ## Run the model and save the result
##' result <- run(model)
##'
##' ## Brief summary of the result. Note that 'U' and 'V' are
##' ## non-empty after running the model.
##' result
setMethod("show",
          signature(object = "SimInf_model"),
          function (object)
          {
              ## The model name
              cat(sprintf("Model: %s\n\n",
                          as.character(class(object))))

              cat(sprintf("Number of nodes: %i\n", Nn(object)))
              cat(sprintf("Number of compartments: %i\n", Nc(object)))
              cat(sprintf("Number of transitions: %i\n", Nt(object)))
              show(object@events)

              cat("\n")
              show_U(object)
              show_V(object)
          }
)

##' Summary of \code{SimInf_model}
##'
##' @aliases summary,SimInf_model-methods
##' @param object The \code{SimInf_model} object
##' @param ... Additional arguments affecting the summary produced.
##' @return None (invisible 'NULL').
##' @keywords methods
##' @export
setMethod("summary",
          signature(object = "SimInf_model"),
          function(object, ...)
          {
              ## The model name
              cat(sprintf("Model: %s\n\n",
                          as.character(class(object))))

              cat(sprintf("Number of nodes: %i\n", Nn(object)))
              cat(sprintf("Number of compartments: %i\n", Nc(object)))
              cat(sprintf("Number of transitions: %i\n", Nt(object)))
              summary(object@events)

              cat("\n")
              show_U(object)
              show_V(object)
          }
)
