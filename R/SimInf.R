## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2026 Stefan Widgren
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

##' A Framework for Data-Driven Stochastic Disease Spread Simulations
##'
##' The \pkg{SimInf} package provides a flexible, high-performance
##' framework for data-driven spatio-temporal disease spread
##' modeling. It is designed to efficiently simulate disease
##' transmission dynamics alongside population demographics and
##' dynamic contact networks.
##'
##' The \pkg{SimInf} framework models infection dynamics within each
##' subpopulation (node) as continuous-time Markov chains (CTMC) using
##' the Gillespie stochastic simulation algorithm (SSA). Additionally,
##' SimInf can incorporate data—such as births, deaths, and
##' movements—as scheduled events. These events trigger at predefined
##' time points and modify the state of subpopulations by randomly
##' sampling individuals from the affected compartments. This
##' capability allows simulations to be driven by empirical records or
##' synthetic scenarios while maintaining the stochastic nature of the
##' population dynamics.
##'
##' The package supports both predefined models (e.g.,
##' \code{\link{SIR}}, \code{\link{SIS}}) and custom model
##' specifications via the \code{\link{mparse}} function.
##' \code{\link{mparse}} serves as the primary interface for defining
##' custom compartment models, allowing users to describe transitions
##' using a simple, human-readable string syntax in R. The function
##' then parses this description, generates model-specific C code, and
##' returns a \code{\linkS4class{SimInf_model}} object ready for
##' simulation. This approach combines the flexibility of R with the
##' computational speed of compiled code, making it well-suited for
##' models with complex propensity functions, multiple compartments,
##' or node-specific parameters.  See the vignette
##' "Getting started with mparse" for a detailed tutorial on defining
##' custom models.
##'
##' After a model is created, a simulation is executed using the
##' \code{\link{run}} method. Upon successful completion,
##' \code{\link{run}} returns a new \code{\linkS4class{SimInf_model}}
##' object containing the original configuration plus the simulated
##' stochastic trajectory.
##'
##' To inspect and analyze the results, SimInf provides a suite of
##' utility functions:
##' \itemize{
##'   \item \code{\link[=summary,SimInf_model-method]{summary}} and
##'     \code{\link[=show,SimInf_model-method]{show}} for a quick
##'     overview of the model structure and simulation results.
##'   \item \code{\link[=plot,SimInf_model-method]{plot}} for
##'     visualizing the time series of compartments and continuous
##'     state variables.
##'   \item \code{\link[=trajectory,SimInf_model-method]{trajectory}}
##'     for extracting the full time series data for custom analysis.
##'   \item \code{\link{prevalence}} for calculating and summarizing
##'     disease prevalence across nodes and time.
##' }
##' These functions facilitate both rapid exploratory analysis and
##' detailed post-processing of simulation outcomes. See the vignette
##' "Post-process data in a trajectory" for a comprehensive tutorial
##' on extracting and analyzing simulation results.
##'
##' Beyond simulation, the package provides functionality to fit
##' models to time series data using two Bayesian inference methods:
##' \itemize{
##'   \item Approximate Bayesian Computation Sequential Monte Carlo
##'     (ABC-SMC), implemented in \code{\link{abc}}, based on the
##'     approach by Toni and others (2009)
##'     \doi{10.1098/rsif.2008.0172}.
##'   \item Particle Markov Chain Monte Carlo (PMCMC), implemented in
##'     \code{\link{pmcmc}}, based on the approach by Andrieu and
##'     others (2010) \doi{10.1111/j.1467-9868.2009.00736.x}.
##' }
##' Both methods enable parameter estimation in stochastic models
##' where the likelihood function is intractable, by using simulated
##' data to estimate the posterior distributions of model parameters.
##'
##' One of our design goal was to make SimInf extendable and enable
##' usage of the numerical solvers from other R extension packages in
##' order to facilitate complex epidemiological research.  To support
##' this, SimInf has functionality to generate the required C and R
##' code from a model specification, see
##' \code{\link{package_skeleton}}
##' @references
##'
##' \Widgren2019
##' @name SimInf
##' @useDynLib SimInf, .registration=TRUE
"_PACKAGE"

##' Unload hook function
##'
##' @param libpath A character string giving the complete path to the
##' package.
##' @noRd
.onUnload <- function(libpath) {
    library.dynam.unload("SimInf", libpath) # nocov
}

##' Example data with spatial distribution of nodes
##'
##' Example data to initialize a population of 1600 nodes and
##' demonstrate various models.
##' @name nodes
##' @docType data
##' @usage data(nodes)
##' @format A \code{data.frame}
##' @keywords dataset
##' @examples
##' \dontrun{
##' ## For reproducibility, call the set.seed() function and specify
##' ## the number of threads to use. To use all available threads,
##' ## remove the set_num_threads() call.
##' set.seed(123)
##' set_num_threads(1)
##'
##' ## Create an 'SIR' model with 1600 nodes and initialize
##' ## it to run over 4*365 days. Add one infected individual
##' ## to the first node.
##' u0 <- u0_SIR()
##' u0$I[1] <- 1
##' tspan <- seq(from = 1, to = 4*365, by = 1)
##' model <- SIR(u0     = u0,
##'              tspan  = tspan,
##'              events = events_SIR(),
##'              beta   = 0.16,
##'              gamma  = 0.077)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##'
##' ## Determine nodes with one or more infected individuals in the
##' ## trajectory. Extract the 'I' compartment and check for any
##' ## infected individuals in each node.
##' infected <- colSums(trajectory(result, ~ I, format = "matrix")) > 0
##'
##' ## Display infected nodes in 'blue' and non-infected nodes in 'yellow'.
##' data("nodes", package = "SimInf")
##' col <- ifelse(infected, "blue", "yellow")
##' plot(y ~ x, nodes, col = col, pch = 20, cex = 2)
##' }
NULL
