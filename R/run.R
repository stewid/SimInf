## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2020 Stefan Widgren
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

##' Keep track of compiled model DLLs.
##' @noRd
.dll <- new.env(parent = emptyenv())

##' Compile the model C code
##'
##' Use 'R CMD SHLIB' to compile the C code for the model and the
##' on-the-fly generated C code to register the native routines for
##' the model.
##' @param model The SimInf model with C code to compile.
##' @param name Character vector with the name of the dll.
##' @param run_fn Name of the function that will be called from
##'     '.Call()'.
##' @return Character vector with the path to the built dll.
##' @noRd
do_compile_model <- function(model, name, run_fn) {
    ## Write the model C code to a temporary file.
    filename <- normalizePath(paste0(tempdir(), "/", name, ".c"),
                              winslash = "/", mustWork = FALSE)
    writeLines(model@C_code, filename)

    ## Include directive for "SimInf.h"
    include <- normalizePath(system.file("include", package = "SimInf"),
                             winslash = "/", mustWork = TRUE)

    ## The output from compiling the model C code.
    lib <- normalizePath(paste0(tempdir(), "/", name, .Platform$dynlib.ext),
                         winslash = "/", mustWork = FALSE)

    ## Set PKG_CPPFLAGS
    pkg_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = NA)
    Sys.setenv(PKG_CPPFLAGS = paste0("-I", shQuote(include),
                                     " -DSIMINF_MODEL_RUN=", run_fn,
                                     " -DSIMINF_R_INIT=R_init_", name,
                                     " -DSIMINF_FORCE_SYMBOLS=FALSE"))

    ## Compile the model C code using the running version of R.
    RBIN <- file.path(R.home(component = "bin"), "R")
    cmd <- paste0(shQuote(RBIN),
                  " CMD SHLIB",
                  " --output=", shQuote(lib),
                  " ", shQuote(filename))
    compiled <- system(cmd, intern = TRUE)

    ## Restore PKG_CPPFLAGS
    if (is.na(pkg_cppflags)) {
        Sys.unsetenv("PKG_CPPFLAGS")
    } else {
        Sys.setenv(PKG_CPPFLAGS = pkg_cppflags)
    }

    if (!file.exists(lib))
        stop(compiled, call. = FALSE)

    lib
}

## Check if model contains C code
contains_C_code <- function(model) {
    if (nchar(paste0(model@C_code, collapse = "\n")))
        return(TRUE)
    FALSE
}

##' Run the SimInf stochastic simulation algorithm
##'
##' @param model The SimInf model to run.
##' @param ... Additional arguments.
##' @param solver Which numerical solver to utilize. Default is 'ssm'.
##' @return \code{\link{SimInf_model}} object with result from simulation.
##' @references \itemize{
##'   \item Bauer P, Engblom S, Widgren S
##'   (2016) "Fast Event-Based Epidemiological Simulations on National Scales"
##'   International Journal of High Performance Computing
##'   Applications, 30(4), 438-453. doi:10.1177/1094342016635723
##'
##'   \item Bauer P., Engblom S. (2015) Sensitivity Estimation and
##'   Inverse Problems in Spatial Stochastic Models of Chemical
##'   Kinetics. In: Abdulle A., Deparis S., Kressner D., Nobile F.,
##'   Picasso M. (eds) Numerical Mathematics and Advanced Applications
##'   - ENUMATH 2013. Lecture Notes in Computational Science and
##'   Engineering, vol 103. Springer, Cham. Doi:
##'   10.1007/978-3-319-10705-9_51
##' }
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
setGeneric("run",
           signature = "model",
           function(model, ...)
               standardGeneric("run"))

##' @rdname run
##' @param solver Which numerical solver to utilize. Default is 'ssm'.
##' @include SimInf_model.R
##' @export
##' @importFrom digest digest
##' @importFrom methods validObject
setMethod("run",
          signature(model = "SimInf_model"),
          function(model, solver = c("ssm", "aem"), ...) {
              solver <- match.arg(solver)
              validObject(model);

              key <- digest(model@C_code, serialize = FALSE)
              if (is.null(.dll[[key]])) {
                  if (!contains_C_code(model))
                      stop("The model must contain C code.")
                  name <- basename(tempfile("SimInf_"))
                  run_fn <- sub("^SimInf_", "run_", name)
                  lib <- do_compile_model(model, name, run_fn)
                  dyn.load(lib)
                  .dll[[key]] <- list(run_fn = run_fn, name = name)
              }

              .Call(.dll[[key]]$run_fn, model, NULL, solver,
                    PACKAGE = .dll[[key]]$name)
          }
)

##' @rdname run
##' @export
setMethod("run",
          signature(model = "SEIR"),
          function(model, solver = c("ssm", "aem"), ...) {
              solver <- match.arg(solver)
              validObject(model);
              .Call(SEIR_run, model, NULL, solver)
          }
)

##' @rdname run
##' @export
setMethod("run",
          signature(model = "SIR"),
          function(model, solver = c("ssm", "aem"), ...) {
              solver <- match.arg(solver)
              validObject(model);
              .Call(SIR_run, model, NULL, solver)
          }
)

##' @rdname run
##' @export
setMethod("run",
          signature(model = "SISe"),
          function(model, solver = c("ssm", "aem"), ...) {
              solver <- match.arg(solver)
              validObject(model);
              .Call(SISe_run, model, NULL, solver)
          }
)

##' @rdname run
##' @export
setMethod("run",
          signature(model = "SISe3"),
          function(model, solver = c("ssm", "aem"), ...) {
              solver <- match.arg(solver)
              validObject(model);
              .Call(SISe3_run, model, NULL, solver)
          }
)

##' @rdname run
##' @export
setMethod("run",
          signature(model = "SISe3_sp"),
          function(model, solver = c("ssm", "aem"), ...) {
              solver <- match.arg(solver)
              validObject(model);
              .Call(SISe3_sp_run, model, NULL, solver)
          }
)

##' @rdname run
##' @export
setMethod("run",
          signature(model = "SISe_sp"),
          function(model, solver = c("ssm", "aem"), ...) {
              solver <- match.arg(solver)
              validObject(model);
              .Call(SISe_sp_run, model, NULL, solver)
          }
)
