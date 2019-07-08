## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015 - 2019  Stefan Engblom
## Copyright (C) 2015 - 2019  Stefan Widgren
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
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Use 'R CMD SHLIB' to compile the C code for the model.
do_compile_model <- function(filename)
{
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

    lib <- paste0(filename, .Platform$dynlib.ext)
    if (!file.exists(lib))
        stop(compiled, call. = FALSE)

    lib
}

## Check if model contains C code
contains_C_code <- function(model)
{
    if (nchar(paste0(model@C_code, collapse = "\n")))
        return(TRUE)
    FALSE
}

##' Run the SimInf stochastic simulation algorithm
##'
##' @param model The siminf model to run.
##' @param threads Number of threads. Default is NULL, i.e. to use all
##'     available processors.
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
           function(model,
                    threads = NULL,
                    solver  = c("ssm", "aem"))
               standardGeneric("run"))

##' @rdname run
##' @include SimInf_model.R
##' @export
##' @importFrom methods validObject
setMethod("run",
          signature(model = "SimInf_model"),
          function(model, threads, solver)
          {
              solver <- match.arg(solver)

              ## FIXME: The 'threads' argument can be dropped with the
              ##  new function 'set_num_threads' added. That also
              ##  means that 'threads' should be removed from the
              ##  expression to parse (below). However, since it is a
              ##  breaking change to remove the 'threads' argument in
              ##  '.Call', just use 'set_num_threads' for now.
              if (!is.null(threads))
                  set_num_threads(threads)

              ## Check that SimInf_model contains all data structures
              ## required by the siminf solver and that they make sense
              validObject(model);

              if (contains_C_code(model)) {
                  ## Write the C code to a temporary file
                  filename <- tempfile("SimInf-")
                  on.exit(unlink(paste0(filename,
                                        c(".c", ".o", .Platform$dynlib.ex))))
                  writeLines(model@C_code, con = paste0(filename, ".c"))

                  lib <- do_compile_model(filename)
                  dll <- dyn.load(lib)
                  on.exit(dyn.unload(lib), add = TRUE)

                  ## Create expression to parse
                  expr <- ".Call(dll$SimInf_model_run, model, threads, solver)"
              } else {
                  ## The model name
                  name <- as.character(class(model))

                  ## The model C run function
                  run_fn <- paste0(name, "_run")

                  ## Create expression to parse
                  expr <- ".Call(run_fn, model, threads, solver, PACKAGE = 'SimInf')"
              }

              ## Run the model. Re-throw any error without the call
              ## included in the error message to make it cleaner.
              tryCatch(eval(parse(text = expr)), error = function(e) {
                  stop(e$message, call. = FALSE)
              })
          }
)
