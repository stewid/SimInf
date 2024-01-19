## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2023 Stefan Widgren
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

## Expression to be evaluated to run a compiled model.
.SimInf_model_run <- paste0(".Call(.dll[[key]]$run_fn, model, solver,",
                            " PACKAGE = .dll[[key]]$name)")

##' Compile the model C code
##'
##' Use 'R CMD SHLIB' to compile the C code for the model and the
##' on-the-fly generated C code to register the native routines for
##' the model.
##' @param model The SimInf model with C code to compile.
##' @param key The digest of the C code to compile.
##' @return Invisible NULL.
##' @noRd
compile_model <- function(model, key) {
    ## Check that the model contains C code.
    if (nchar(paste0(model@C_code, collapse = "\n")) == 0)
        stop("The model must contain C code.", call. = FALSE)

    ## Determine the name and run_fun to call from R.
    name <- basename(tempfile("SimInf_"))
    run_fn <- sub("^SimInf_", "run_", name)

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
    compiled <- tools::Rcmd(c("SHLIB",
                              paste0("--output=", shQuote(lib)),
                              shQuote(filename)),
                            stdout = TRUE,
                            stderr = TRUE)

    ## Restore PKG_CPPFLAGS
    if (is.na(pkg_cppflags)) {
        Sys.unsetenv("PKG_CPPFLAGS")
    } else {
        Sys.setenv(PKG_CPPFLAGS = pkg_cppflags)
    }

    if (!file.exists(lib))
        stop(compiled, call. = FALSE)

    dyn.load(lib)
    .dll[[key]] <- list(run_fn = run_fn, name = name)

    invisible(NULL)
}

##' Determine the key to the compiled model DLL.
##' @noRd
model_dll_key <- function(model) {
    key <- digest::digest(model@C_code)
    if (is.null(.dll[[key]]))
        compile_model(model, key)
    key
}

##' Run the SimInf stochastic simulation algorithm
##'
##' @param model The SimInf model to run.
##' @param ... Additional arguments.
##' @param solver Which numerical solver to utilize. Default is 'ssm'.
##' @return \code{\link{SimInf_model}} object with result from
##'     simulation.
##' @references
##'
##' \Widgren2019
##'
##' \Bauer2016
##'
##' \Bauer2015
##' @examples
##' ## For reproducibility, call the set.seed() function and specify
##' ## the number of threads to use. To use all available threads,
##' ## remove the set_num_threads() call.
##' set.seed(123)
##' set_num_threads(1)
##'
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
setGeneric(
    "run",
    signature = "model",
    function(model, ...) {
        standardGeneric("run")
    }
)

##' @rdname run
##' @param solver Which numerical solver to utilize. Default is 'ssm'.
##' @include SimInf_model.R
##' @export
setMethod(
    "run",
    signature(model = "SimInf_model"),
    function(model, solver = c("ssm", "aem"), ...) {
        solver <- match.arg(solver)
        methods::validObject(model)
        key <- model_dll_key(model)
        eval(parse(text = .SimInf_model_run))
    }
)

##' @rdname run
##' @export
setMethod(
    "run",
    signature(model = "SEIR"),
    function(model, solver = c("ssm", "aem"), ...) {
        solver <- match.arg(solver)
        methods::validObject(model)
        .Call(SEIR_run, model, solver)
    }
)

##' @rdname run
##' @export
setMethod(
    "run",
    signature(model = "SIR"),
    function(model, solver = c("ssm", "aem"), ...) {
        solver <- match.arg(solver)
        methods::validObject(model)
        .Call(SIR_run, model, solver)
    }
)

##' @rdname run
##' @export
setMethod(
    "run",
    signature(model = "SIS"),
    function(model, solver = c("ssm", "aem"), ...) {
        solver <- match.arg(solver)
        methods::validObject(model)
        .Call(SIS_run, model, solver)
    }
)

##' @rdname run
##' @export
setMethod(
    "run",
    signature(model = "SISe"),
    function(model, solver = c("ssm", "aem"), ...) {
        solver <- match.arg(solver)
        methods::validObject(model)
        .Call(SISe_run, model, solver)
    }
)

##' @rdname run
##' @export
setMethod(
    "run",
    signature(model = "SISe3"),
    function(model, solver = c("ssm", "aem"), ...) {
        solver <- match.arg(solver)
        methods::validObject(model)
        .Call(SISe3_run, model, solver)
    }
)

##' @rdname run
##' @export
setMethod(
    "run",
    signature(model = "SISe3_sp"),
    function(model, solver = c("ssm", "aem"), ...) {
        solver <- match.arg(solver)
        methods::validObject(model)
        .Call(SISe3_sp_run, model, solver)
    }
)

##' @rdname run
##' @export
setMethod(
    "run",
    signature(model = "SISe_sp"),
    function(model, solver = c("ssm", "aem"), ...) {
        solver <- match.arg(solver)
        methods::validObject(model)
        .Call(SISe_sp_run, model, solver)
    }
)
