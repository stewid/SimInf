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

##' Function to compile custom \code{\link{SimInf_model}}.
##'
##' This function compiles a model specified using \code{\link{mparse}}
##' function, and produces a customised named \code{\link{SimInf_model}}
##' object that can be run in the usual way. This is useful for routines
##' that require multiple calls to the \code{run} method for 
##' \code{\link{SimInf_model}} objects created using \code{\link{mparse}},
##' which avoids the need to re-compile the model each time \code{run} is
##' called.
##' @include SimInf_model.R
##' @export
##' @examples
##' ## Create an SIR model object using \code{mparse}.
##' transitions <- c("S -> beta*S*I -> I", "I -> gamma*I -> R")
##' compartments <- c("S", "I", "R")
##' u0 <- data.frame(S = 99, I = 1, R = 0)
##' model <- mparse(transitions = transitions, compartments = compartments,
##'     gdata = c(beta = 0.16, gamma = 0.077), u0 = u0, tspan = 1:100)
##'
##' ## Run the SIR model and plot the result.
##' set.seed(22)
##' result <- run(model)
##' plot(result)
##' 
##' ## compile the model first and then re-run
##' ## and plot the result
##' set.seed(22)
##' model <- compile_mparse(model, "SIR")
##' result <- run(model)
##' plot(result)
compile_mparse <- function(model, name) {

    ## Check that SimInf_model contains all data structures
    ## required by the siminf solver and that they make sense
    validObject(model)
    
    if(missing(name)) {
        stop("No 'name' argument provided")
    }
    if(!is.character(name)) {
        stop("'name' is not character")
    }
    if(length(name) != 1) {
        stop("'name' not character of length 1")
    }

    if (nchar(paste0(model@C_code, collapse = "\n"))) {
        ## Write the C code to a temporary file
        filename <- paste0("SimInf-", name)
        unlink(paste0(filename,
                            c(".c", ".o", .Platform$dynlib.ex)))
        writeLines(model@C_code, con = paste0(filename, ".c"))

        ## Include directive for "SimInf.h"
        include <- system.file("include", package = "SimInf")
        Sys.setenv(PKG_CPPFLAGS=sprintf("-I%s", shQuote(include)))

        ## Compile the model C code using the running version of R.
        cmd <- paste(shQuote(file.path(R.home(component = "bin"), "R")),
                   "CMD SHLIB",
                   shQuote(paste0(basename(filename), ".c")))
        compiled <- system(cmd, intern = TRUE)

        ## check compilation
        lib <- paste0(filename, .Platform$dynlib.ext)
        if (!file.exists(lib)) {
            stop(compiled)
        }
                
        ## update C_code slot
        model@C_code <- character()
    } else {
        stop("No C code to compile")
    }
    ## output new class and model
    setClass(filename, contains = c("SimInf_model"), where = topenv(parent.frame()))
    as(model, filename)
}

