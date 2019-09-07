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

do_parse_prior <- function(prior) {
    prior <- as.character(prior)
    if (!identical(length(prior), 3L)) {
        stop("Invalid formula specification for prior.",
             call. = FALSE)
    }

    ## Determine the parameter to fit.
    parameter <- prior[2]

    ## Determine the distribution for the parameter.
    pattern <- paste0("^([GNU])\\s*\\(\\s*",
                      "([-+]?[0-9]*\\.?[0-9]+),\\s*",
                      "([-+]?[0-9]*\\.?[0-9]+)\\)$")
    d <- prior[3]
    m <- regexec(pattern, d)
    if (m[[1]][1] == -1) {
        stop("Invalid formula specification for priors.",
             call. = FALSE)
    }
    m <- regmatches(d, m)[[1]][-1]
    distribution <- m[1]
    hyperparameters <- as.numeric(m[-1])

    ## Check hyperparameters.
    if (distribution == "G" && !all(hyperparameters > 0)) {
        stop("Invalid prior: gamma hyperparameters must be > 0.",
             call. = FALSE)
    }
    if (distribution == "N" && hyperparameters[2] < 0) {
        stop("Invalid prior: normal variance must be > 0.",
             call. = FALSE)
    }
    if (distribution == "U" && hyperparameters[1] >= hyperparameters[2]) {
        stop("Invalid prior: uniform bounds in wrong order.",
             call. = FALSE)
    }

    data.frame(parameter = parameter, distribution = distribution,
               p1 = hyperparameters[1], p2 = hyperparameters[2],
               stringsAsFactors = FALSE)
}

##' @noRd
##' @importFrom stats dgamma
##' @importFrom stats dnorm
##' @importFrom stats dunif
##' @importFrom stats rgamma
##' @importFrom stats rnorm
##' @importFrom stats runif
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
    priors <- do.call("rbind", lapply(priors, do_parse_prior))

    if (any(duplicated(priors$parameter)) || any(nchar(priors$parameter) == 0))
        stop("'priors' must have non-duplicated parameter names.")

    priors
}
