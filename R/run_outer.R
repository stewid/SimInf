## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
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

##' Run \code{SimInf_model} on scaled parameters
##'
##' @param x Scale the model \code{gdata} parameter values on the
##'     right hand side of the formula with \code{x} before calling
##'     \code{FUN} with the scaled model as argument.
##' @param y Scale the model \code{gdata} parameter values on the
##'     left hand side of the formula with \code{y} before calling
##'     \code{FUN} with the scaled model as argument.
##' @param model The siminf model to scale parameters on and run.
##' @param formula The parameters in the \code{gdata} vector matching
##'     the left hand side of the formula \code{a + b ~ c} will be
##'     scaled by \code{y}.  The parameters in the \code{gdata} vector
##'     matching the right hand side of the formula \code{a + b ~ c}
##'     will be scaled by \code{x}.
##' @param FUN A function to use on the scaled model 'gdata' parameters.
##' @param ... Optional arguments to be passed to \code{FUN}.
##' @return Array with dimension \code{c(dim(x), dim(y))}.
##' @include SimInf_model.R
##' @export
##' @importFrom stats terms
##' @examples
##' ## Create an SIR-model with 100 nodes of 99 susceptible individuals
##' ## and one infected individuals.
##' u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
##' model <- SIR(u0, 1:75, beta = 0.16, gamma = 0.077)
##'
##' ## Define scaling parameters
##' x <- seq(from = 0.2, to = 1.8, by = 0.1)
##' y <- seq(from = 0.2, to = 1.1, by = 0.1)
##'
##' ## Utility function to run the model and estimate the population
##' ## prevalence on day 75.
##' pop_prev <- function(model) {
##'     result <- run(model)
##'     prevalence(result, I~., type = "pop", as.is = TRUE)[75]
##' }
##'
##' ## Scale 'gamma' with 'y' and 'beta' with 'x' and
##' ## run the model and determine the population prevalence on day
##' ## 500. For each combination of 'x' and 'y', the model parameters
##' ## are scaled and the function 'pop_prev' called with the
##' ## perturbed model.
##' pop <- run_outer(x, y, model, gamma ~ beta, pop_prev)
##'
##' ## Plot result
##' contour(x * model@gdata["beta"], y * model@gdata["gamma"],
##'         pop, method = "edge", bty = "l")
##'
##' ## Utility function to run the model and estimate the node
##' ## prevalence on day 75.
##' node_prev <- function(model) {
##'     result <- run(model)
##'     prevalence(result, I~., type = "nop", as.is = TRUE)[75]
##' }
##'
##' ## Scale 'gamma' with 'y' and 'beta' with 'x' and
##' ## run the model and determine the node prevalence on day
##' ## 500. For each combination of 'x' and 'y', the model parameters
##' ## are scaled and the function 'node_prev' called with the
##' ## perturbed model.
##' nop <- run_outer(x, y, model, gamma ~ beta, node_prev)
##'
##' ## Plot result
##' contour(x * model@gdata["beta"], y * model@gdata["gamma"],
##'         nop, method = "edge", bty = "l")
run_outer <- function(x, y, model, formula = NULL, FUN = NULL, ...)
{
    ## Check 'x'
    if (missing(x))
        stop("Missing 'x' argument.", call. = FALSE)
    if (!is.numeric(x))
        stop("'x' argument is not numeric.", call. = FALSE)

    ## Check 'y'
    if (missing(y))
        stop("Missing 'y' argument.", call. = FALSE)
    if (!is.numeric(y))
        stop("'y' argument is not numeric.", call. = FALSE)

    check_model_argument(model)

    if (is.null(names(model@gdata)))
        stop("'names(model@gdata)' is NULL.", call. = FALSE)
    if (is.null(formula))
        stop("'formula' argument is NULL.", call. = FALSE)
    if (is.null(FUN))
        stop("'FUN' argument is NULL.", call. = FALSE)
    FUN <- match.fun(FUN)

    ## Determine indices to the 'gdata' parameters to scale by 'x'
    xx <- attr(terms(formula, allowDotAsName = TRUE), "term.labels")
    xx <- xx[attr(terms(formula, allowDotAsName = TRUE), "order") == 1]
    if (length(xx) < 1)
        stop("Invalid parameters on the right side of the formula.", call. = FALSE)
    x_i <- match(xx, names(model@gdata))
    if (any(is.na(x_i)))
        stop("Unmatched parameters on the right side of the formula.", call. = FALSE)

    ## Determine indices to the 'gdata' parameters to scale by 'y'
    yy <- attr(terms(formula, allowDotAsName = TRUE), "response")
    if (yy < 1)
        stop("Invalid parameters on the left side of the formula.", call. = FALSE)
    vars <- attr(terms(formula, allowDotAsName = TRUE), "variables")[-1]
    yy <- as.character(vars[yy])
    yy <- unlist(strsplit(yy, "+", fixed = TRUE))
    yy <- sub("^\\s", "", sub("\\s$", "", yy))
    y_i <- match(yy, names(model@gdata))
    if (any(is.na(y_i))) {
        stop("Unmatched parameters on the left hand side of the formula.",
             call. = FALSE)
    }

    outer(x, y, function(x, y, ...) {
        run_internal <- function(x, y, x_i, y_i, model, ...) {
            model@gdata[x_i] <- model@gdata[x_i] * x
            model@gdata[y_i] <- model@gdata[y_i] * y
            FUN(model, ...)
        }

        sapply(seq_len(length(x)), function(i, ...) {
            run_internal(x[i], y[i], x_i, y_i, model, ...)
        }, ...)
    }, ...)
}
