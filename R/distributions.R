## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
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

check_binomial_distribution <- function(hyperparameters, symbols) {
    if (length(hyperparameters) != 2 || anyNA(hyperparameters)) {
        stop("Invalid formula specification for binomial distribution.",
             call. = FALSE)
    }

    if (is.null(symbols)) {
        if (!all(is.numeric(hyperparameters), all(hyperparameters >= 0))) {
            stop("Invalid distribution: binomial hyperparameters must be >= 0.",
                 call. = FALSE)
        }

        if (!is_wholenumber(hyperparameters[1])) {
            stop("Invalid distribution: ",
                 "binomial size must be an integer >= 0.",
                 call. = FALSE)
        }

        if (hyperparameters[2] > 1) {
            stop("Invalid distribution: ",
                 "binomial probability must be <= 1.",
                 call. = FALSE)
        }
    }
}

check_gamma_distribution <- function(hyperparameters, symbols) {
    if (length(hyperparameters) != 2 || anyNA(hyperparameters)) {
        stop("Invalid formula specification for gamma distribution.",
             call. = FALSE)
    }

    if (is.null(symbols) && !all(hyperparameters > 0)) {
        stop("Invalid distribution: ",
             "gamma hyperparameters must be > 0.",
             call. = FALSE)
    }
}

check_normal_distribution <- function(hyperparameters, symbols) {
    if (length(hyperparameters) != 2 || anyNA(hyperparameters)) {
        stop("Invalid formula specification for normal distribution.",
             call. = FALSE)
    }

    if (is.null(symbols) && hyperparameters[2] < 0) {
        stop("Invalid distribution: ",
             "normal variance must be >= 0.",
             call. = FALSE)
    }
}

check_poisson_distribution <- function(hyperparameters, symbols) {
    if (length(hyperparameters) != 2 ||
        is.na(hyperparameters[1]) ||
        !is.na(hyperparameters[2])) {
            stop("Invalid formula specification for poisson distribution.",
             call. = FALSE)
    }

    if (is.null(symbols) && hyperparameters[1] < 0) {
        stop("Invalid distribution: ",
             "lambda must be >= 0.",
             call. = FALSE)
    }
}

check_uniform_distribution <- function(hyperparameters, symbols) {
    if (length(hyperparameters) != 2 || anyNA(hyperparameters)) {
        stop("Invalid formula specification for uniform distribution.",
             call. = FALSE)
    }

    if (is.null(symbols) && hyperparameters[1] >= hyperparameters[2]) {
        stop("Invalid distribution: ",
             "uniform bounds in wrong order.",
             call. = FALSE)
    }
}

check_hyperparameters <- function(distribution, hyperparameters, symbols) {
    switch(distribution,
           binomial = {
               check_binomial_distribution(hyperparameters, symbols)
           },
           gamma = {
               check_gamma_distribution(hyperparameters, symbols)
           },
           normal = {
               check_normal_distribution(hyperparameters, symbols)
           },
           poisson = {
               check_poisson_distribution(hyperparameters, symbols)
           },
           uniform = {
               check_uniform_distribution(hyperparameters, symbols)
           },
           stop("Unknown distribution: '", distribution, "'.", call. = FALSE)
           )
}

##' Determine if there are any symbols in the hyperparameter(s).
##' @noRd
parse_hyperparameter_symbols <- function(tokens) {
    unlist(lapply(seq_len(nrow(tokens)), function(i) {
        if (tokens$token[i] == "SYMBOL")
            return(tokens$text[i])
        NULL
    }))
}

##' Determine the hyperparameter(s) for the distribution.
##' @noRd
parse_hyperparameters <- function(distribution, tokens, symbols) {
    comma <- which(tokens$token == "','")
    if (length(comma) == 0) {
        hyperparameters <- c(paste0(tokens$text, collapse = ""), NA)
    } else if (length(comma) == 1) {
        hyperparameters <- c(paste0(tokens$text[seq_len(comma - 1)],
                                    collapse = ""),
                             paste0(tokens$text[-seq_len(comma)],
                                    collapse = ""))
    } else {
        stop("Invalid formula specification for distribution.",
             call. = FALSE)
    }

    if (is.null(symbols))
        hyperparameters <- suppressWarnings(as.numeric(hyperparameters))

    check_hyperparameters(distribution, hyperparameters, symbols)

    hyperparameters
}

parse_distribution <- function(dist) {
    err_str <- "Invalid formula specification for distribution."
    dist <- as.character(dist)
    if (!identical(length(dist), 3L))
        stop(err_str, call. = FALSE)

    ## Determine the parameter from the lhs.
    parameter <- dist[2]

    ## Parse the rhs of the formula.
    tokens <- utils::getParseData(parse(text = dist[3], keep.source = TRUE))
    tokens <- tokens[tokens$terminal, c("token", "text")]
    if (!all(tokens$token[1] == "SYMBOL_FUNCTION_CALL",
             tokens$token[2] == "'('",
             tokens$token[nrow(tokens)] == "')'")) {
        stop(err_str, call. = FALSE)
    }

    ## Determine the distribution for the parameter.
    distribution <- tokens$text[1]
    tokens <- tokens[c(-1, -2, -nrow(tokens)), 1:2]

    ## Determine the hyperparameter(s) for the distribution.
    symbols <- parse_hyperparameter_symbols(tokens)
    hyperparameters <- parse_hyperparameters(distribution, tokens, symbols)

    data.frame(parameter = parameter, distribution = distribution,
               p1 = hyperparameters[1], p2 = hyperparameters[2],
               symbols = I(list(symbols)), stringsAsFactors = FALSE)
}

##' @noRd
parse_priors <- function(priors) {
    if (is.list(priors)) {
        if (!all(vapply(priors, is, logical(1), "formula"))) {
            stop("'priors' must be a formula or a list with formula items.",
                 call. = FALSE)
        }
    } else if (!methods::is(priors, "formula")) {
        stop("'priors' must be a formula or a list with formula items.",
             call. = FALSE)
    } else {
        priors <- list(priors)
    }

    ## Determine priors for parameters in the model
    priors <- do.call("rbind", lapply(priors, function(prior) {
        prior <- parse_distribution(prior)
        if (!is.null(unlist(prior$symbols))) {
            stop("Invalid formula specification for 'priors'.",
                 call. = FALSE)
        }
        prior[, c("parameter", "distribution", "p1", "p2")]
    }))

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
    if (anyNA(pars)) {
        pars <- match(priors$parameter, names(model@gdata))
        if (anyNA(pars)) {
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

##' Generate random deviates from priors.
##' @noRd
rpriors <- function(priors, n = 1) {
    mapply(function(parameter, distribution, p1, p2) {
        switch(
            distribution,
            "gamma" = stats::rgamma(n = n, shape = p1, rate = 1 / p2),
            "normal" = stats::rnorm(n = n, mean = p1, sd = p2),
            "uniform" = stats::runif(n = n, min = p1, max = p2),
            stop("Unknown distribution.", call. = FALSE))
    },
    priors$parameter,
    priors$distribution,
    priors$p1,
    priors$p2)
}

##' Determine the sum of the log of the densities for the priors.
##' @noRd
dpriors <- function(x, priors) {
    sum(mapply(function(x, distribution, p1, p2) {
        switch(
            distribution,
            "gamma" = stats::dgamma(x = x,
                                    shape = p1,
                                    rate = 1 / p2,
                                    log = TRUE),
            "normal" = stats::dnorm(x = x,
                                    mean = p1,
                                    sd = p2,
                                    log = TRUE),
            "uniform" = stats::dunif(x = x,
                                     min = p1,
                                     max = p2,
                                     log = TRUE),
            stop("Unknown distribution.", call. = FALSE))
    },
    x,
    priors$distribution,
    priors$p1,
    priors$p2))
}
