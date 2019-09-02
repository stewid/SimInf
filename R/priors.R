## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2019 Stefan Widgren
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

##' @noRd
parse_priors <- function(priors) {
    if (is.list(priors)) {
        if (!all(vapply(priors, is, logical(1), "formula"))) {
            stop("'priors' must be a formula or a list with formula items.",
                 call. = FALSE)
        }
    } else if (!is(priors, "formula")) {
        stop("'priors' must be a formula or a list with formula items.",
             call. = FALSE)
    } else {
        priors <- list(priors)
    }

    ## Determine priors for parameters in the model
    parameters <- lapply(priors, function(f) {
        f <- as.character(f)
        if (!identical(length(f), 3L)) {
            stop("Invalid formula specification for prior.",
                 call. = FALSE)
        }

        ## Determine the parameter to fit.
        parameter <- f[2]

        ## Determine the distribution for the parameter.
        pattern <- paste0("^([GNU])\\s*\\(\\s*",
                          "([-+]?[0-9]*\\.?[0-9]+),\\s*",
                          "([-+]?[0-9]*\\.?[0-9]+)\\)$")
        d <- f[3]
        m <- regexec(pattern, d)
        if (m[[1]][1] == -1)
            stop("Invalid formula specification for priors.", call. = FALSE)
        m <- regmatches(d, m)[[1]][-1]
        distribution <- m[1]
        hyperparameters <- as.numeric(m[-1])

        ## Check hyperparameters.
        if (distribution == "G") {
            if (!all(hyperparameters > 0)) {
                stop("Invalid prior: gamma hyperparameters must be > 0.",
                     call. = FALSE)
            }
        } else if (distribution == "N") {
            if (hyperparameters[2] < 0) {
                stop("Invalid prior: normal variance must be > 0.",
                     call. = FALSE)
            }
        } else if (distribution == "U") {
            if (hyperparameters[1] >= hyperparameters[2]) {
                stop("Invalid prior: uniform bounds in wrong order.",
                     call. = FALSE)
            }
        } else {
            stop("Invalid prior: unknown distribution.", call. = FALSE)
        }

        ## Generate random numbers according to the distribution.
        rfun <- switch(distribution, G = rgamma, N = rnorm, U = runif)
        random <- function(n) {
            rfun(n, hyperparameters[1], hyperparameters[2])
        }

        ## Gives the density for the distribution.
        dfun <- switch(distribution, G = dgamma, N = dnorm, U = dunif)
        density <- function(x, log = FALSE) {
            dfun(x, hyperparameters[1], hyperparameters[2], log)
        }

        list(parameter       = parameter,
             distribution    = distribution,
             hyperparameters = hyperparameters,
             random          = random,
             density         = density)
    })

    lbl <- vapply(parameters, function(x) x$parameter, character(1))
    if (is.null(lbl) || any(duplicated(lbl)) || any(nchar(lbl) == 0))
        stop("'priors' must have non-duplicated parameter names.")
    names(parameters) <- lbl
    parameters
}
