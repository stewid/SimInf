## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015 - 2018  Stefan Engblom
## Copyright (C) 2015 - 2018  Stefan Widgren
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

    ## Check 'type' argument
    type <- match.arg(type)

    ## Check the 'node' argument
    if (!is.null(node)) {
        if (!is.numeric(node))
            stop("'node' must be integer")
        if (!all(is_wholenumber(node)))
            stop("'node' must be integer")
        if (min(node) < 1)
            stop("'node' must be integer > 0")
        if (max(node) > Nn(model))
            stop("'node' must be integer <= number of nodes")
        node <- as.integer(sort(unique(node)))
    }

    ## Determine compartments for population
    j <- attr(terms(formula, allowDotAsName = TRUE), "term.labels")
    j <- j[attr(terms(formula, allowDotAsName = TRUE), "order") == 1]
    if (length(j) < 1)
        stop("Invalid formula specification of population")
    pop <- unlist(sapply(j, function(jj) {
        ## Replace '.' with all discrete compartments in the model.
        if (identical(jj, "."))
            jj <- rownames(model@S)
        jj
    }))
    pop <- unique(as.character(pop))
    if (!length(pop))
        stop("'pop' is empty")
    j <- !(pop %in% rownames(model@S))
    if (any(j)) {
        stop("Non-existing compartment(s) in model: ",
             paste0("'", pop[j], "'", collapse = ", "))
    }

    ## Determine compartments for cases
    j <- attr(terms(formula, allowDotAsName = TRUE), "response")
    if (j < 1)
        stop("Invalid formula specification of 'cases'")
    cases <- attr(terms(formula, allowDotAsName = TRUE), "variables")[-1]
    j <- as.character(cases[j])
    j <- unlist(strsplit(j, "+", fixed = TRUE))
    j <- sub("^\\s", "", sub("\\s$", "", j))
    cases <- unlist(sapply(j, function(jj) {
        ## Replace '.' with all discrete compartments in the model.
        if (identical(jj, "."))
            jj <- rownames(model@S)
        jj
    }))
    cases <- unique(as.character(cases))
    if (!length(cases))
        stop("'cases' is empty")
    j <- !(cases %in% rownames(model@S))
    if (any(j)) {
        stop("Non-existing compartment(s) in model: ",
             paste0("'", cases[j], "'", collapse = ", "))
    }

    ## Sum all individuals in 'cases' compartments in a matrix with
    ## one row per node X length(tspan)
    cm <- NULL
    for (compartment in cases) {
        if (is.null(cm)) {
            cm <- trajectory(model, compartments = compartment,
                             node = node, as.is = TRUE)
        } else {
            cm <- cm + trajectory(model, compartments = compartment,
                                  node = node, as.is = TRUE)
        }
    }
    dimnames(cm) <- NULL

    ## Sum all individuals in 'pop' compartments in a matrix with one
    ## row per node X length(tspan)
    pm <- NULL
    for (compartment in pop) {
        if (is.null(pm)) {
            pm <- trajectory(model, compartments = compartment,
                             node = node, as.is = TRUE)
        } else {
            pm <- pm + trajectory(model, compartments = compartment,
                                  node = node, as.is = TRUE)
        }
    }
    dimnames(pm) <- NULL

    if (identical(type, "pop")) {
        cm <- colSums(cm)
        pm <- colSums(pm)
    } else if (identical(type, "nop")) {
        cm <- colSums(cm > 0)
        ## Only include nodes with individuals
        pm <- colSums(pm > 0)
    }

    if (isTRUE(as.is))
        return(cm / pm)

    time <- names(model@tspan)
    if (is.null(time))
        time <- model@tspan
    if (type %in% c("pop", "nop")) {
        return(data.frame(time = time,
                          prevalence = cm / pm,
                          stringsAsFactors = FALSE))
    }

    if (is.null(node))
        node = seq_len(Nn(model))

    data.frame(node = node,
               time = rep(time, each = length(node)),
               prevalence = as.numeric(cm / pm),
               stringsAsFactors = FALSE)
}