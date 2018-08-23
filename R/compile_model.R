## SimInf, a framework for stochastic disease spread simulations
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

##' Class to handle a compiled custom \code{\link{SimInf_model}}
##'
##' @section Slots:
##' \describe{
##'   \item{filename}{
##'     Character vector of length 1 containing name of shared library.
##'   }
##' }
##' @include SimInf_model.R
##' @export
setClass("SimInf_model_dll",
         contains = c("SimInf_model"),
         slots = c(filename = "character"),
         validity = function(object) {
             ## check filename
             if (!is.character(object@filename) | length(object@filename) != 1)
                return("'filename' is not a character of length 1")

             TRUE
         }
)

##' @rdname run
##' @value \code{\link{SimInf_model_dll}} with result from simulation.
##' @export
setMethod("run",
          signature(model = "SimInf_model_dll"),
          function(model, threads, solver)
          {
              solver <- match.arg(solver)

              ## Check that SimInf_model contains all data structures
              ## required by the siminf solver and that they make sense
              validObject(model);

              dll <- dyn.load(model@filename)
              on.exit(dyn.unload(model@filename), add = TRUE)

              ## Create expression to parse
              expr <- ".Call(dll$SimInf_model_run, model, threads, solver)"

              ## Run the model. Re-throw any error without the call
              ## included in the error message to make it cleaner.
              tryCatch(eval(parse(text = expr)), error = function(e) {
                  stop(e$message, call. = FALSE)
              })
          }
)

##' Function to compile custom \code{\link{SimInf_model}}.
##'
##' This function compiles a model specified using \code{\link{mparse}}
##' function, and produces a \code{\link{SimInf_model}} object with an
##' internal link to a shared library file that can be run in the usual way.
##' This is useful for routines that require multiple calls to the \code{run}
##' method for \code{\link{SimInf_model}} objects, since it avoids the need
##' to re-compile the model each time \code{run} is called.
##' @param model     An object of class \code{\link{SimInf_model}}.
##' @param filename A character specifying the name of the shared
##'     library that will be created. Default is NULL, i.e. to use a
##'     temporary file.
##' @include SimInf_model.R
##' @return \linkS4class{SimInf_model_dll}
##' @export
##' @importFrom methods as
##' @examples
##' ## Create an SIR model object using mparse.
##' transitions <- c("S -> beta*S*I -> I", "I -> gamma*I -> R")
##' compartments <- c("S", "I", "R")
##' u0 <- data.frame(S = 99, I = 1, R = 0)
##' model <- mparse(transitions = transitions,
##'     compartments = compartments,
##'     gdata = c(beta = 0.16, gamma = 0.077),
##'     u0 = u0, tspan = 1:100)
##'
##' ## Run the SIR model and plot the result.
##' ## This recompiles the model when run()
##' ## is called
##' set.seed(22)
##' result <- run(model)
##' plot(result)
##'
##' ## Compile the model first and then re-run
##' set.seed(22)
##' model <- compile_model(model, "SIR")
##' class(model)
##' result <- run(model)
##' plot(result)
compile_model <- function(model, filename = NULL)
{
    check_model_argument(model)

    if (is.null(filename))
        filename <- tempfile("SimInf-")
    if (!is.character(filename) || length(filename) != 1)
        stop("'filename' not character of length 1")

    if (!contains_C_code(model))
        stop("No C code to compile")

    ## Write the C code to a temporary file
    unlink(paste0(filename, c(".c", ".o", .Platform$dynlib.ex)))
    writeLines(model@C_code, con = paste0(filename, ".c"))

    ## output revised model
    model <- as(model, "SimInf_model_dll")
    model@filename <- do_compile_model(filename)
    model
}
