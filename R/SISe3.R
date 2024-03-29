## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
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

##' Definition of the \sQuote{SISe3} model
##'
##' Class to handle the SISe3 \code{\link{SimInf_model}} model.
##' @include SimInf_model.R
##' @export
setClass("SISe3", contains = c("SimInf_model"))

##' The compartments in an SISe3 model
##' @noRd
compartments_SISe3 <- function() {
    c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3")
}

##' The select matrix 'E' for an SISe3 model
##' @noRd
select_matrix_SISe3 <- function() {
    matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
             1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1),
           nrow = 6,
           ncol = 6,
           dimnames = list(compartments_SISe3(), seq_len(6)))
}

##' Create a \code{SISe3} model
##'
##' Create a \code{SISe3} model to be used by the simulation
##' framework.
##'
##' The \code{SISe3} model contains two compartments in three age
##' categories; number of susceptible (S_1, S_2, S_3) and number of
##' infectious (I_1, I_2, I_3). Additionally, it contains an
##' environmental compartment to model shedding of a pathogen to the
##' environment. Consequently, the model has six state transitions,
##'
##' \deqn{S_1 \stackrel{\upsilon_1 \varphi S_1}{\longrightarrow} I_1}{
##' S_1 -- upsilon_1 phi S_1 --> I_1}
##'
##' \deqn{I_1 \stackrel{\gamma_1 I_1}{\longrightarrow} S_1}{
##' I_1 -- gamma_1 I_1 --> S_1}
##'
##' \deqn{S_2 \stackrel{\upsilon_2 \varphi S_2}{\longrightarrow} I_2}{
##' S_2 -- upsilon_2 phi S_2 --> I_2}
##'
##' \deqn{I_2 \stackrel{\gamma_2 I_2}{\longrightarrow} S_2}{
##' I_2 -- gamma_2 I_2 --> S_2}
##'
##' \deqn{S_3 \stackrel{\upsilon_3 \varphi S_3}{\longrightarrow} I_3}{
##' S_3 -- upsilon_3 phi S_3 --> I_3}
##'
##' \deqn{I_3 \stackrel{\gamma_3 I_3}{\longrightarrow} S_3}{
##' I_3 -- gamma_3 I_3 --> S_3}
##'
##' where the transition rate per unit of time from susceptible to
##' infected is proportional to the concentration of the environmental
##' contamination \eqn{\varphi}{phi} in each node. Moreover, the
##' transition rate from infected to susceptible is the recovery rate
##' \eqn{\gamma_1, \gamma_2, \gamma_3}, measured per individual and
##' per unit of time. Finally, the environmental infectious pressure
##' in each node is evolved by,
##'
##' \deqn{\frac{d\varphi(t)}{dt} = \frac{\alpha \left(I_1(t) + I_2(t)
##' + I_3(t)\right)}{N(t)} - \beta(t) \varphi(t) + \epsilon}{
##' dphi(t) / dt = alpha (I_1(t) + I_2(t) + I_3(t)) / N(t)
##' - beta(t) phi(t) + epsilon}
##'
##' where \eqn{\alpha} is the average shedding rate of the pathogen to
##' the environment per infected individual and \eqn{N = S_1 + S_2 +
##' S_3 + I_1 + I_2 + I_3} the size of the node. The seasonal decay
##' and removal of the pathogen is captured by \eqn{\beta(t)}. It is
##' also possible to include a small background infectious pressure
##' \eqn{\epsilon} to allow for other indirect sources of
##' environmental contamination. The environmental infectious pressure
##' \eqn{\varphi(t)}{phi(t)} in each node is evolved each time unit by
##' the Euler forward method. The value of \eqn{\varphi(t)}{phi(t)} is
##' saved at the time-points specified in \code{tspan}.
##'
##' The argument \code{u0} must be a \code{data.frame} with one row for
##' each node with the following columns:
##' \describe{
##' \item{S_1}{The number of sucsceptible in age category 1}
##' \item{I_1}{The number of infected in age category 1}
##' \item{S_2}{The number of sucsceptible in age category 2}
##' \item{I_2}{The number of infected in age category 2}
##' \item{S_3}{The number of sucsceptible in age category 3}
##' \item{I_3}{The number of infected in age category 3}
##' }
##'
##' @template beta-section
##' @template u0-param
##' @template tspan-param
##' @template events-param
##' @template phi-param
##' @param upsilon_1 Indirect transmission rate of the environmental
##' infectious pressure in age category 1
##' @param upsilon_2 Indirect transmission rate of the environmental
##' infectious pressure in age category 2
##' @param upsilon_3 Indirect transmission rate of the environmental
##' infectious pressure in age category 3
##' @param gamma_1 The recovery rate from infected to susceptible for
##' age category 1
##' @param gamma_2 The recovery rate from infected to susceptible for
##' age category 2
##' @param gamma_3 The recovery rate from infected to susceptible for
##' age category 3
##' @param alpha Shed rate from infected individuals
##' @template beta-end-param
##' @param epsilon The background environmental infectious pressure
##' @return \code{SISe3}
##' @include check_arguments.R
##' @export
SISe3 <- function(u0,
                  tspan,
                  events    = NULL,
                  phi       = NULL,
                  upsilon_1 = NULL,
                  upsilon_2 = NULL,
                  upsilon_3 = NULL,
                  gamma_1   = NULL,
                  gamma_2   = NULL,
                  gamma_3   = NULL,
                  alpha     = NULL,
                  beta_t1   = NULL,
                  beta_t2   = NULL,
                  beta_t3   = NULL,
                  beta_t4   = NULL,
                  end_t1    = NULL,
                  end_t2    = NULL,
                  end_t3    = NULL,
                  end_t4    = NULL,
                  epsilon   = NULL) {

    ## Check arguments.

    ## Check u0 and compartments
    u0 <- check_u0(u0, compartments_SISe3())

    ## Check initial infectious pressure
    if (is.null(phi))
        phi <- 0
    phi <- rep(phi, length.out = nrow(u0))
    check_infectious_pressure_arg(nrow(u0), phi)

    ## Check 'gdata' parameters
    check_gdata_arg(upsilon_1, upsilon_2, upsilon_3, gamma_1, gamma_2, gamma_3,
                    alpha, beta_t1, beta_t2, beta_t3, beta_t4, epsilon)

    ## Check interval endpoints
    check_integer_arg(end_t1, end_t2, end_t3, end_t4)
    end_t1 <- rep(end_t1, length.out = nrow(u0))
    end_t2 <- rep(end_t2, length.out = nrow(u0))
    end_t3 <- rep(end_t3, length.out = nrow(u0))
    end_t4 <- rep(end_t4, length.out = nrow(u0))
    check_end_t_arg(nrow(u0), end_t1, end_t2, end_t3, end_t4)

    ## Arguments seem ok...go on

    N <- matrix(c(2, 2, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, -1, 0, -1, 0, -1),
                nrow = 6, ncol = 3,
                dimnames = list(compartments_SISe3(), c("1", "2", "3")))

    G <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                  0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1),
                nrow   = 6, ncol   = 6,
                dimnames = list(c("S_1 -> upsilon_1*phi*S_1 -> I_1",
                                  "I_1 -> gamma_1*I_1 -> S_1",
                                  "S_2 -> upsilon_2*phi*S_2 -> I_2",
                                  "I_2 -> gamma_2*I_2 -> S_2",
                                  "S_3 -> upsilon_3*phi*S_3 -> I_3",
                                  "I_3 -> gamma_3*I_3 -> S_3"),
                                c("1", "2", "3", "4", "5", "6")))

    S <- matrix(c(-1, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0,
                  0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, -1),
                nrow = 6, ncol = 6,
                dimnames = list(compartments_SISe3(),
                                c("1", "2", "3", "4", "5", "6")))

    v0 <- matrix(as.numeric(phi), nrow  = 1, byrow = TRUE,
                 dimnames = list("phi"))

    ldata <- matrix(as.numeric(c(end_t1, end_t2, end_t3, end_t4)),
                    nrow = 4, byrow = TRUE,
                    dimnames = list(c("end_t1", "end_t2", "end_t3", "end_t4")))

    gdata <- as.numeric(c(upsilon_1, upsilon_2, upsilon_3,
                          gamma_1, gamma_2, gamma_3, alpha,
                          beta_t1, beta_t2, beta_t3, beta_t4, epsilon))
    names(gdata) <- c("upsilon_1", "upsilon_2", "upsilon_3",
                      "gamma_1", "gamma_2", "gamma_3", "alpha",
                      "beta_t1", "beta_t2", "beta_t3", "beta_t4", "epsilon")

    model <- SimInf_model(G      = G,
                          S      = S,
                          E      = select_matrix_SISe3(),
                          N      = N,
                          tspan  = tspan,
                          events = events,
                          ldata  = ldata,
                          gdata  = gdata,
                          u0     = u0,
                          v0     = v0)

    methods::as(model, "SISe3")
}

##' Example data to initialize events for the \sQuote{SISe3} model
##'
##' Example data to initialize scheduled events for a population of
##' 1600 nodes and demonstrate the \code{\linkS4class{SISe3}} model.
##'
##' Example data to initialize scheduled events (see
##' \code{\linkS4class{SimInf_events}}) for a population of 1600 nodes
##' and demonstrate the \code{\linkS4class{SISe3}} model. The dataset
##' contains 783773 events for 1600 nodes distributed over 4 * 365
##' days. The events are divided into three types: \sQuote{Exit}
##' events remove individuals from the population (n = 182535),
##' \sQuote{Enter} events add individuals to the population (n =
##' 182685), \sQuote{Internal transfer} events move individuals
##' between compartmens within one node e.g. ageing (n = 317081), and
##' \sQuote{External transfer} events move individuals between nodes
##' in the population (n = 101472). The vignette contains a detailed
##' description of how scheduled events operate on a model.
##' @name events_SISe3
##' @docType data
##' @usage data(events_SISe3)
##' @format A \code{data.frame}
##' @keywords dataset
##' @examples
##' ## For reproducibility, call the set.seed() function and specify
##' ## the number of threads to use. To use all available threads,
##' ## remove the set_num_threads() call.
##' set.seed(123)
##' set_num_threads(1)
##'
##' ## Create an 'SISe3' model with 1600 nodes and initialize
##' ## it to run over 4*365 days. Add one infected individual
##' ## to the first node.
##' data("u0_SISe3", package = "SimInf")
##' data("events_SISe3", package = "SimInf")
##' u0_SISe3$I_1[1] <- 1
##' tspan <- seq(from = 1, to = 4*365, by = 1)
##' model <- SISe3(u0 = u0_SISe3, tspan = tspan, events = events_SISe3,
##'                phi = rep(0, nrow(u0_SISe3)), upsilon_1 = 1.8e-2,
##'                upsilon_2 = 1.8e-2, upsilon_3 = 1.8e-2,
##'                gamma_1 = 0.1, gamma_2 = 0.1, gamma_3 = 0.1,
##'                alpha = 1, beta_t1 = 1.0e-1, beta_t2 = 1.0e-1,
##'                beta_t3 = 1.25e-1, beta_t4 = 1.25e-1, end_t1 = 91,
##'                end_t2 = 182, end_t3 = 273, end_t4 = 365, epsilon = 0)
##'
##' ## Display the number of individuals affected by each event type
##' ## per day.
##' plot(events(model))
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##'
##' ## Summarize the trajectory. The summary includes the number of
##' ## events by event type.
##' summary(result)
NULL

##' Example data to initialize the \sQuote{SISe3} model
##'
##' Example data to initialize a population of 1600 nodes and
##' demonstrate the \code{\linkS4class{SISe3}} model.
##'
##' A \code{data.frame} with the number of individuals in the
##' \sQuote{S_1}, \sQuote{S_2}, \sQuote{S_3}, \sQuote{I_1},
##' \sQuote{I_2} and \sQuote{I_3} compartments in 1600 nodes. Note
##' that the \sQuote{I_1}, \sQuote{I_2} and \sQuote{I_3} compartments
##' are zero.
##' @name u0_SISe3
##' @docType data
##' @usage data(u0_SISe3)
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
##' ## Create an 'SISe3' model with 1600 nodes and initialize it to
##' ## run over 4*365 days and record data at weekly time-points.
##'
##' ## Load the initial population and add ten infected individuals to
##' ## I_1 in the first node.
##' u0 <- u0_SISe3
##' u0$I_1[1] <- 10
##'
##' ## Define 'tspan' to run the simulation over 4*365 and record the
##' ## state of the system at weekly time-points.
##' tspan <- seq(from = 1, to = 4*365, by = 7)
##'
##' ## Load scheduled events for the population of nodes with births,
##' ## deaths and between-node movements of individuals.
##' events <- events_SISe3
##'
##' ## Create a 'SISe3' model
##' model <- SISe3(u0 = u0, tspan = tspan, events = events,
##'                phi = rep(0, nrow(u0)), upsilon_1 = 1.8e-2,
##'                upsilon_2 = 1.8e-2, upsilon_3 = 1.8e-2,
##'                gamma_1 = 0.1, gamma_2 = 0.1, gamma_3 = 0.1,
##'                alpha = 1, beta_t1 = 1.0e-1, beta_t2 = 1.0e-1,
##'                beta_t3 = 1.25e-1, beta_t4 = 1.25e-1, end_t1 = 91,
##'                end_t2 = 182, end_t3 = 273, end_t4 = 365, epsilon = 0)
##'
##' ## Run the model to generate a single stochastic trajectory.
##' result <- run(model)
##'
##' ## Summarize trajectory
##' summary(result)
##'
##' ## Plot the proportion of nodes with at least one infected
##' ## individual.
##' plot(result, I_1 + I_2 + I_3 ~ ., level = 2, type = "l")
##' }
NULL
