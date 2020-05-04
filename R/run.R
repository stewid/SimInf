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

##' Compile the model C code
##'
##' Use 'R CMD SHLIB' to compile the C code for the model and the
##' on-the-fly generated C code to register the native routines for
##' the model.
##' @param model The SimInf model with C code to compile.
##' @param name Character vector with the name of the dll.
##' @return Character vector with the path to the built dll.
##' @noRd
do_compile_model <- function(model, name) {
    lines <- c(
        "#include <Rdefines.h>",
        "#include <R_ext/Rdynload.h>",
        "#include <R_ext/Visibility.h>",
        "",
        "SEXP SimInf_model_run(SEXP, SEXP, SEXP);",
        "",
        "static const R_CallMethodDef callMethods[] =",
        "{",
        "    {\"SimInf_model_run\", (DL_FUNC)&SimInf_model_run, 3},",
        "    {NULL, NULL, 0}",
        "};",
        "",
        paste0("void attribute_visible R_init_", name, "(DllInfo *info)"),
        "{",
        "    R_registerRoutines(info, NULL, callMethods, NULL, NULL);",
        "    R_useDynamicSymbols(info, TRUE);",
        "    R_forceSymbols(info, FALSE);",
        "}",
        "")

    ## Write the model C code to a temporary file.
    filename <- file.path(tempdir(), paste0(name, ".c"))
    writeLines(model@C_code, filename)

    ## Write the model init C code to a temporary file.
    filename_init <- file.path(tempdir(), paste0(name, "_init.c"))
    writeLines(lines, filename_init)

    ## Include directive for "SimInf.h"
    include <- system.file("include", package = "SimInf")
    Sys.setenv(PKG_CPPFLAGS = sprintf("-I%s", shQuote(include)))

    ## Compile the model C code using the running version of R.
    wd <- setwd(tempdir())
    cmd <- paste(shQuote(file.path(R.home(component = "bin"), "R")),
                 "CMD SHLIB",
                 shQuote(basename(filename)),
                 shQuote(basename(filename_init)))
    compiled <- system(cmd, intern = TRUE)
    setwd(wd)

    lib <- file.path(tempdir(), paste0(name, .Platform$dynlib.ext))
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

              ## Check that SimInf_model contains all data structures
              ## required by the siminf solver and that they make sense
              validObject(model);

              if (contains_C_code(model)) {
                  name <- paste0("SimInf_",
                                 digest(model@C_code, serialize = FALSE))
                  if (!is.loaded("SimInf_model_run", name, "Call")) {
                      lib <- do_compile_model(model, name)
                      dyn.load(lib)
                  }

                  ## Run the model
                  return(.Call("SimInf_model_run", model, NULL,
                               solver, PACKAGE = name))
              }

              ## The model name
              name <- as.character(class(model))

              ## The model C run function
              run_fn <- paste0(name, "_run")

              ## Run the model
              eval(parse(text = ".Call(run_fn, model, NULL, solver)"))
          }
)
