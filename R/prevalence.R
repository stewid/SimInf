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

## Sum all individuals in compartments in a matrix with one row per
## node X length(tspan)
sum_compartments <- function(model, compartments, index) {
    m <- NULL

    for (j in seq_len(length(compartments))) {
        for (compartment in names(compartments[[j]])) {
            if (is.null(m)) {
                m <- trajectory(model, compartment, index, "matrix")
            } else {
                m <- m + trajectory(model, compartment, index, "matrix")
            }
        }
    }

    m
}

evaluate_condition <- function(model, compartments, index, n) {
    ## Create an environment to hold the trajectory data with one
    ## column for each compartment.
    e <- new.env(parent = baseenv())
    for (j in seq_len(length(compartments$rhs))) {
        if (length(compartments$rhs[[j]]) > 0) {
            ac <- attr(compartments$rhs[[j]], "available_compartments")
            for (compartment in ac) {
                assign(x = compartment,
                       value = as.integer(
                           trajectory(model, compartment, index, "matrix")),
                       pos = e)
            }
        }
    }

    ## Then evaluate the condition using the data in the environment.
    condition <- compartments$condition
    e$condition <- condition
    k <- evalq(eval(parse(text = condition)), envir = e)
    l <- length(model@tspan) * ifelse(is.null(index), n, length(index))
    if (!is.logical(k) || length(k) != l) {
        stop(paste0("The condition must be either 'TRUE' ",
                    "or 'FALSE' for every node and time step."),
             call. = FALSE)
    }

    matrix(k, ncol = length(model@tspan))
}

calculate_prevalence <- function(model, compartments, level,
                                 index, n, format, id) {
    ## Sum all individuals in the 'cases' and 'population'
    ## compartments in a matrix with one row per node X length(tspan)
    cases <- sum_compartments(model, compartments$lhs, index)
    population <- sum_compartments(model, compartments$rhs, index)

    ## Apply condition
    if (!is.null(compartments$condition)) {
        condition <- evaluate_condition(model, compartments, index, n)
        cases <- cases * condition
        population <- population * condition
    }

    if (identical(level, 1L)) {
        cases <- colSums(cases)
        population <- colSums(population)
    } else if (identical(level, 2L)) {
        cases <- colSums(cases > 0)
        ## Only include nodes with individuals
        population <- colSums(population > 0)
    }

    prevalence <- cases / population

    if (identical(format, "matrix"))
        return(prevalence)

    time <- names(model@tspan)
    if (is.null(time))
        time <- model@tspan
    if (level %in% c(1L, 2L)) {
        return(data.frame(time = time,
                          prevalence = prevalence,
                          stringsAsFactors = FALSE))
    }

    if (is.null(index))
        index <- seq_len(n)

    prevalence <- data.frame(id = index,
                             time = rep(time, each = length(index)),
                             prevalence = as.numeric(prevalence),
                             stringsAsFactors = FALSE)
    colnames(prevalence)[1] <- id
    prevalence
}

##' Calculate prevalence from a model object with trajectory data
##'
##' Calculate the proportion of individuals with disease in the
##' population, or the proportion of nodes with at least one diseased
##' individual, or the proportion of individuals with disease in each
##' node.
##' @param model The \code{model} with trajectory data to calculate
##'     the prevalence from.
setGeneric(
    "prevalence",
    signature = c("model", "formula"),
    function(model,
             formula,
             level = 1,
             index = NULL,
             format = c("data.frame", "matrix")) {
        standardGeneric("prevalence")
    }
)

##' @rdname prevalence
##' @param formula A formula that specifies the compartments that
##'     define the cases with a disease or that have a specific
##'     characteristic (numerator), and the compartments that define
##'     the entire population of interest (denominator). The
##'     left-hand-side of the formula defines the cases, and the
##'     right-hand-side defines the population, for example,
##'     \code{I~S+I+R} in a \sQuote{SIR} model (see
##'     \sQuote{Examples}). The \code{.}  (dot) is expanded to all
##'     compartments, for example, \code{I~.}  is expanded to
##'     \code{I~S+I+R} in a \sQuote{SIR} model (see
##'     \sQuote{Examples}). The formula can also contain a condition
##'     (indicated by \code{|}) for each node and time step to further
##'     control the population to include in the calculation, for
##'     example, \code{I ~ . | R == 0} to calculate the prevalence
##'     when the recovered is zero in a \sQuote{SIR} model. The
##'     condition must evaluate to \code{TRUE} or \code{FALSE} in each
##'     node and time step. Note that if the denominator is zero, the
##'     prevalence is \code{NaN}.
##' @param level The level at which the prevalence is calculated at
##'     each time point in \code{tspan}. 1 (population prevalence):
##'     calculates the proportion of the individuals (cases) in the
##'     population. 2 (node prevalence): calculates the proportion of
##'     nodes with at least one case. 3 (within-node prevalence):
##'     calculates the proportion of cases within each node. Default
##'     is \code{1}.
##' @param index Indices specifying the subset of nodes to include in
##'     the calculation of the prevalence. Default is \code{index =
##'     NULL}, which includes all nodes.
##' @param format The default (\code{format = "data.frame"}) is to
##'     generate a \code{data.frame} with one row per time-step with
##'     the prevalence. Using \code{format = "matrix"} returns the
##'     result as a matrix.
##' @return A \code{data.frame} if \code{format = "data.frame"}, else
##'     a matrix.
##' @include SimInf_model.R
##' @include match_compartments.R
##' @export
##' @examples
##' ## Create an 'SIR' model with 6 nodes and initialize
##' ## it to run over 10 days.
##' u0 <- data.frame(S = 100:105, I = c(0, 1, 0, 2, 0, 3), R = rep(0, 6))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
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
##' prevalence(result, I~S+I+R, level = 2)
##'
##' ## Determine the proportion of infected individuals in each node
##' ## at the time-points in 'tspan'.
##' prevalence(result, I~S+I+R, level = 3)
##'
##' ## Determine the proportion of infected individuals in each node
##' ## at the time-points in 'tspan' when the number of recovered is
##' ## zero.
##' prevalence(result, I~S+I+R|R==0, level = 3)
setMethod(
    "prevalence",
    signature(model = "SimInf_model", formula = "formula"),
    function(model, formula, level, index, format) {
        compartments <- match_compartments(compartments = formula,
                                           ok_combine = FALSE,
                                           ok_lhs = TRUE,
                                           U = rownames(model@S))
        if (is.null(compartments$lhs))
            stop("Invalid 'formula' specification.", call. = FALSE)

        check_integer_arg(level)
        level <- as.integer(level)
        if (length(level) != 1 || any(level < 1) || any(level > 3))
            stop("'level' must be an integer with a value 1, 2 or 3.")
        index <- check_node_index_argument(model, index)
        format <- match.arg(format)
        n <- n_nodes(model)
        id <- "node"

        calculate_prevalence(model, compartments, level, index, n, format, id)
    }
)
