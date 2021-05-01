## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2021 Stefan Widgren
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

##' @importFrom utils getParseData
##' @noRd
do_parse_prior <- function(prior) {
    err_str <- "Invalid formula specification for distribution."
    prior <- as.character(prior)
    if (!identical(length(prior), 3L))
        stop(err_str, call. = FALSE)

    ## Determine the parameter to fit from the lhs.
    parameter <- prior[2]

    ## Parse the rhs of the formula.
    tokens <- getParseData(parse(text = prior[3], keep.source = TRUE))
    tokens <- tokens[tokens$terminal, c("token", "text")]
    if (!all(tokens$token[1] == "SYMBOL_FUNCTION_CALL",
             tokens$token[2] == "'('",
             tokens$token[nrow(tokens)] == "')'")) {
        stop(err_str, call. = FALSE)
    }

    ## Determine the distribution for the parameter.
    distribution <- tokens$text[1]
    tokens <- tokens[c(-1, -2, -nrow(tokens)), 1:2]

    ## Determine the hyperparameters for the distribution.
    comma <- which(tokens$token == "','")
    if (length(comma) != 1)
        stop(err_str, call. = FALSE)
    hyperparameters <- c(paste0(tokens$text[seq_len(comma - 1)], collapse = ""),
                         paste0(tokens$text[-seq_len(comma)], collapse = ""))
    hyperparameters <- suppressWarnings(as.numeric(hyperparameters))
    if (any(length(hyperparameters) != 2, any(is.na(hyperparameters))))
        stop(err_str, call. = FALSE)

    ## Check distribution and hyperparameters.
    switch(distribution,
           gamma = {
               if (!all(hyperparameters > 0)) {
                   stop("Invalid distribution: ",
                        "gamma hyperparameters must be > 0.",
                        call. = FALSE)
               }
           },
           normal = {
               if (hyperparameters[2] < 0) {
                   stop("Invalid distribution: ",
                        "normal variance must be > 0.",
                        call. = FALSE)
               }
           },
           uniform = {
               if (hyperparameters[1] >= hyperparameters[2]) {
                   stop("Invalid distribution: ",
                        "uniform bounds in wrong order.",
                        call. = FALSE)
               }
           },
           stop("'distribution' must be one of 'gamma', 'normal' or 'uniform'.",
                call. = FALSE)
           )

    data.frame(parameter = parameter, distribution = distribution,
               p1 = hyperparameters[1], p2 = hyperparameters[2],
               stringsAsFactors = FALSE)
}

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
    priors <- do.call("rbind", lapply(priors, do_parse_prior))

    if (any(duplicated(priors$parameter)) ||
        any(nchar(priors$parameter) == 0)) {
        stop("'priors' must have non-duplicated parameter names.",
             call. = FALSE)
    }

    priors
}

##' Match the 'priors' to parameters in 'ldata' or 'gdata'.
##' @noRd
match_priors <- function(model, priors) {
    pars <- match(priors$parameter, rownames(model@ldata))
    if (any(is.na(pars))) {
        pars <- match(priors$parameter, names(model@gdata))
        if (any(is.na(pars))) {
            stop("All parameters in 'priors' must be either ",
                 "in 'gdata' or 'ldata'.", call. = FALSE)
        }
        target <- "gdata"
    } else {
        if (!identical(n_nodes(model), 1L))
            stop("The 'model' must contain one node.", call. = FALSE)
        target <- "ldata"
    }

    list(pars = pars, target = target)
}
