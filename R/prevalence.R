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

## Determine the compartments in the formula item and split
## 'compartment1 + compartment2 + ...'. Moreover, trim whitespace and
## replace '.' with all compartments in the model.
parse_formula_item <- function(x, compartments)
{
    x <- unlist(strsplit(x, "+", fixed = TRUE))
    x <- sub("^\\s", "", sub("\\s$", "", x))
    x <- unlist(sapply(x, function(y) {
        if (identical(y, "."))
            y <- compartments
        y
    }))
    x <- unique(as.character(x))
    if (!length(x))
        stop("No compartments in formula specification.")
    i <- !(x %in% compartments)
    if (any(i)) {
        stop("Non-existing compartment(s) in model: ",
             paste0("'", x[i], "'", collapse = ", "))
    }
    x
}

## Sum all individuals in compartments in a matrix with one row per
## node X length(tspan)
sum_individuals <- function(model, compartments, node)
{
    m <- NULL
    for (compartment in compartments) {
        if (is.null(m)) {
            m <- trajectory(model, compartments = compartment,
                            node = node, as.is = TRUE)
        } else {
            m <- m + trajectory(model, compartments = compartment,
                                node = node, as.is = TRUE)
        }
    }
    dimnames(m) <- NULL
    m
}

##' Calculate prevalence from a model object with trajectory data
##'
##' Calculate the proportion of individuals with disease in the
##' population, or the proportion of nodes with at least one diseased
##' individual, or the proportion of individuals with disease in each
##' node.
##' @param model The \code{model} with trajectory data to calculate
##'     the prevalence from.
##' @param formula A formula that specify the compartments that define
##'     the cases with a disease or a condition (numerator), and the
##'     compartments that define the entire population of interest
##'     (denominator). The left hand side of the formula defines the
##'     cases, and the right hand side defines the population, for
##'     example, \code{I~S+I+R} in a \sQuote{SIR} model (see
##'     \sQuote{Examples}). The \code{.}  (dot) is expanded to all
##'     compartments, for example, \code{I~.}  is expanded to
##'     \code{I~S+I+R} in a \sQuote{SIR} model (see
##'     \sQuote{Examples}).
##' @param type The type of prevalence measure to calculate at each
##'     time point in \code{tspan}: \code{pop} (population prevalence)
##'     calculates the proportion of the individuals (cases) in the
##'     population, \code{nop} (node prevalence) calculates the
##'     proportion of nodes with at least one case, and \code{wnp}
##'     (within-node prevalence) calculates the proportion of cases
##'     within each node. Default is \code{pop}.
##' @param node Indices specifying the subset nodes to include in the
##'     calculation of the prevalence. Default is \code{NULL}, which
##'     includes all nodes.
##' @param as.is The default (\code{as.is = FALSE}) is to generate a
##'     \code{data.frame} with one row per time-step with the
##'     prevalence. Using \code{as.is = TRUE} returns the result as a
##'     matrix, which is the internal format.
##' @return A \code{data.frame} if \code{as.is = FALSE}, else a
##'     matrix.
##' @include SimInf_model.R
##' @include check_arguments.R
##' @export
##' @importFrom stats terms
##' @examples
##' ## Create an 'SIR' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = c(0, 1, 0, 2, 0, 3), R = rep(0, 6))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model, threads = 1)
##'
##' ## Determine the proportion of infected individuals (cases)
##' ## in the population at the time-points in 'tspan'.
##' prevalence(result, I~S+I+R)
##'
##' ## Identical result is obtained with the shorthand 'I~.'
##' prevalence(result, I~.)
##'
##' ## Determine the proportion of nodes with infected individuals at
##' ## the time-points in 'tspan'.
##' prevalence(result, I~S+I+R, type = "nop")
##'
##' ## Determine the proportion of infected individuals in each node
##' ## at the time-points in 'tspan'.
##' prevalence(result, I~S+I+R, type = "wnp")
prevalence <- function(model,
                       formula,
                       type = c("pop", "nop", "wnp"),
                       node = NULL,
                       as.is = FALSE)
{
    check_model_argument(model)

    ## Check 'formula' argument
    if (missing(formula))
        stop("Missing 'formula' argument")
    if (!is(formula, "formula"))
        stop("'formula' argument is not a 'formula'")
    formula <- as.character(formula)
    if (!identical(length(formula), 3L))
        stop("Invalid formula specification.")

    ## Check 'type' argument
    type <- match.arg(type)

    ## Check the 'node' argument
    node <- check_node_argument(model, node)

    ## Determine the compartments for the cases and the population.
    cases <- parse_formula_item(formula[2], rownames(model@S))
    population <- parse_formula_item(formula[3], rownames(model@S))

    ## Sum all individuals in the 'cases' and 'population'
    ## compartments in a matrix with one row per node X length(tspan)
    cases <- sum_individuals(model, cases, node)
    population <- sum_individuals(model, population, node)

    if (identical(type, "pop")) {
        cases <- colSums(cases)
        population <- colSums(population)
    } else if (identical(type, "nop")) {
        cases <- colSums(cases > 0)
        ## Only include nodes with individuals
        population <- colSums(population > 0)
    }

    prevalence <- cases / population

    if (isTRUE(as.is))
        return(prevalence)

    time <- names(model@tspan)
    if (is.null(time))
        time <- model@tspan
    if (type %in% c("pop", "nop")) {
        return(data.frame(time = time,
                          prevalence = prevalence,
                          stringsAsFactors = FALSE))
    }

    if (is.null(node))
        node = seq_len(Nn(model))

    data.frame(node = node,
               time = rep(time, each = length(node)),
               prevalence = as.numeric(prevalence),
               stringsAsFactors = FALSE)
}
