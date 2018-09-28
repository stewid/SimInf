## SimInf, a framework for stochastic disease spread simulations
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

##' Class \code{"SimInf_controller"}
##'
##' Class to handle a SimInf controller. The goal of the controller is
##' to make it possible to interact and probe data when running a
##' trajectory with a solver.
##' @section Slots: \describe{ \item{model}{ The
##'     \code{linkS4class{SimInf_model}} object to apply the
##'     controller on.  } \item{C_code}{ Character vector with
##'     optional controller C code. If non-empty, the C code is
##'     written to a temporary C-file when the \code{run} method is
##'     called.  The temporary C-file is compiled and the resulting
##'     DLL is dynamically loaded. The DLL is unloaded and the
##'     temporary files are removed after running the controller.  } }
##' @export
##' @importFrom methods validObject
setClass("SimInf_controller",
         slots = c(model  = "SimInf_model",
                   C_code = "character"),
         validity = function(object) {
             ## Check model
             errors <- validObject(object@model)
             if (!isTRUE(errors))
                 return(errors)

             TRUE
         }
)

##' Create a \code{SimInf_controller}
##'
##' @param model The \code{linkS4class{SimInf_model}} object to apply
##'     the controller on.
##' @param C_code Character vector with optional controller C code. If
##'     non-empty, the C code is written to a temporary C-file when
##'     the \code{run} method is called.  The temporary C-file is
##'     compiled and the resulting DLL is dynamically loaded. The DLL
##'     is unloaded and the temporary files are removed after running
##'     the controller.
##' @return \linkS4class{SimInf_controller} object.
##' @export
##' @importFrom methods new
SimInf_controller <- function(model  = NULL, C_code = NULL)
{
    new("SimInf_controller", model = model, C_code = C_code)
}

##' @rdname run
##' @export
##' @importFrom methods validObject
setMethod("run",
          signature(model = "SimInf_controller"),
          function(model, threads, solver)
          {
              solver <- match.arg(solver)

              ## Check that model contains all required data
              ## structures and that they make sense.
              validObject(model);

              if (!contains_C_code(model@model))
                  stop("The 'model' object must contain C code")
              if (!contains_C_code(model))
                  stop("The 'controller' object must contain C code")

              ## Write the model C code to a temporary file
              filename <- tempfile("SimInf-")
              on.exit(unlink(paste0(filename, ".c")), add = TRUE)
              on.exit(unlink(paste0(filename, ".o")), add = TRUE)
              on.exit(unlink(paste0(filename, .Platform$dynlib.ex)), add = TRUE)
              writeLines(model@model@C_code, con = paste0(filename, ".c"))

              lib_model <- do_compile_model(filename)
              dll_model <- dyn.load(lib_model)
              on.exit(dyn.unload(lib_model), add = TRUE)

              ## Write the controller C code to a temporary file
              filename <- tempfile("SimInf-")
              on.exit(unlink(paste0(filename, ".c")), add = TRUE)
              on.exit(unlink(paste0(filename, ".o")), add = TRUE)
              on.exit(unlink(paste0(filename, .Platform$dynlib.ex)), add = TRUE)
              writeLines(model@C_code, con = paste0(filename, ".c"))

              lib_controller <- do_compile_model(filename)
              dll_controller <- dyn.load(lib_controller)
              on.exit(dyn.unload(lib_controller), add = TRUE)

              ## Create expression to parse
              expr <- paste0(".Call(",
                             "dll_controller$SimInf_controller_run, ",
                             "dll_model$SimInf_model_run$address, ",
                             "model@model, ",
                             "threads, ",
                             "solver)")

              ## Run controller and model.
              eval(parse(text = expr))
          }
)
