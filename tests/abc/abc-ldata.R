## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2022 Stefan Widgren
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

library(SimInf)
library(tools)
source("util/check.R")

## Specify the number of threads to use.
set_num_threads(1)

## Create a model with parameters in ldata
model <- mparse(transitions = c("S -> beta*S*I/(S+I+R) -> I + Icum",
                                "I -> gamma*I -> R"),
                compartments = c("S", "I", "Icum", "R"),
                ldata = data.frame(beta = 1, gamma = 0.5),
                u0 = data.frame(S = 9999, I = 0, Icum = 0, R = 0),
                events = data.frame(event = 1, time = 25, node = 1,
                                    dest = 0, n = 1, proportion = 0,
                                    select = 1, shift = 0),
                E = matrix(c(0, 1, 0, 0), nrow = 4, ncol = 1,
                           dimnames = list(c("S", "I", "Icum", "R"),
                                           c("1"))),
                tspan = 2:75)

## Check that a non-numeric distance vector raises an error.
res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = 2,
                       distance = function(result, ...) {
                           c("1", "2")
                       },
                       tolerance = c(5, 4)))
check_error(res, "The result from the ABC distance function must be numeric.")

## Check that a distance vector with the wrong dimension raises an
## error.
res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = 2,
                       distance = function(result, ...) {
                           1:3
                       },
                       tolerance = c(5, 4)))
check_error(
    res,
    "Invalid dimension of the result from the ABC distance function.")

## Check that a distance vector with NA raises an error.
res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = 2,
                       distance = function(result, ...) {
                           c(NA, 2:20)
                       },
                       tolerance = c(5, 4)))
check_error(
    res,
    "The result from the ABC distance function must be non-negative.")

## Check that a distance vector with a negative value raises an error.
res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = 2,
                       distance = function(result, ...) {
                           c(-1L, 2:20)
                       },
                       tolerance = c(5, 4)))
check_error(
    res,
    "The result from the ABC distance function must be non-negative.")

## Check that a non-numeric tolerance raises an error.
res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = 2,
                       distance = function(result, ...) {
                           1:20
                       },
                       tolerance = c("1", "2")))
check_error(res, "'tolerance' must have non-negative values.")

## Check that a tolerance with NA raises an error.
res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = 2,
                       distance = function(result, ...) {
                           1:20
                       },
                       tolerance = c(NA_real_, 2)))
check_error(res, "'tolerance' must have non-negative values.")

## Check that a 'tolerance' with a negative value raises an error.
res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = 2,
                       distance = function(result, ...) {
                           1:20
                       },
                       tolerance = c(1, -2)))
check_error(res, "'tolerance' must have non-negative values.")

## Check that a tolerance with wrong dimension raises an error.
res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = 2,
                       distance = function(result, ...) {
                           1:20
                       },
                       tolerance = numeric(0)))
check_error(res, "'tolerance' must have columns.")

## Check that a non-decreasing tolerance raises an error.
res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = 2,
                       distance = function(result, ...) {
                           1:20
                       },
                       tolerance = c(4, 5)))
check_error(res, "'tolerance' must be a decreasing vector.")

distance_fn_ldata <- function(result, ...) {
    ## Extract the time-series for R1 for each node as a data.frame.
    sim <- trajectory(result, "Icum")

    ## Split the 'sim' data.frame by node and calculate the sum of the
    ## squared distance at each time-point for every node.
    tapply(sim$Icum, sim$node, function(Icum) {
        ## Observed cases
        cases <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 6, 6, 25, 42, 56,
                   106, 171, 279, 382, 576, 710, 977, 934, 846, 672,
                   585, 430, 346, 221, 192, 172, 122, 66, 48, 57, 26,
                   12, 10, 6, 6, 8, 5, 0, 1, 2, 1, 0, 4, 0, 0, 0, 1,
                   0, 0, 0, 0, 0, 0)

        ## Simulated cases
        sim_cases <- c(0, diff(c(0, Icum)))

        sum((sim_cases - cases)^2)
    })
}

## Check invalid npart
res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = 1,
                       distance = distance_fn_ldata,
                       tolerance = c(250000, 225000)))
check_error(res, "'npart' must be an integer > 1.")

res <- assertError(abc(model = model,
                       priors = c(beta ~ uniform(0.5, 1.5),
                                  gamma ~ uniform(0.3, 0.7)),
                       npart = c(10, 10),
                       distance = distance_fn_ldata,
                       tolerance = c(250000, 225000)))
check_error(res, "'npart' must be an integer > 1.")

set.seed(123)
fit <- abc(model = model,
           priors = c(beta ~ uniform(0.5, 1.5),
                      gamma ~ uniform(0.3, 0.7)),
           npart = 10,
           distance = distance_fn_ldata,
           tolerance = c(250000, 225000),
           verbose = TRUE)
fit
summary(fit)
as.data.frame(fit)

fit <- continue(fit, tolerance = 200000, verbose = TRUE)

pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(fit, xlim = c(0.3, 1.5), ylim = c(0.3, 1.5))
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that an invalid 'n' is detected.
sigma <- SimInf:::abc_proposal_covariance(SimInf:::abc_particles(fit, 2L))
res <- assertError(
    .Call(SimInf:::SimInf_abc_proposals,
          fit@priors$parameter,
          fit@priors$distribution,
          fit@priors$p1,
          fit@priors$p2,
          0L,
          SimInf:::abc_particles(fit, 2L),
          fit@weight[, 2],
          sigma))
check_error(res, "'n' must be an integer > 0.")

res <- assertError(
    .Call(SimInf:::SimInf_abc_proposals,
          fit@priors$parameter,
          fit@priors$distribution,
          fit@priors$p1,
          fit@priors$p2,
          2,
          SimInf:::abc_particles(fit, 2L),
          fit@weight[, 2],
          sigma))
check_error(res, "'n' must be an integer > 0.")

## Check that an invalid 'parameter' is detected.
res <- assertError(
    .Call(SimInf:::SimInf_abc_proposals,
          3L,
          fit@priors$distribution,
          fit@priors$p1,
          fit@priors$p2,
          1L,
          SimInf:::abc_particles(fit, 2L),
          fit@weight[, 2],
          sigma))
check_error(res, "'parameter' must be a character vector.")

## Check that an invalid 'distribution' is detected.
res <- assertError(
    .Call(SimInf:::SimInf_abc_proposals,
          fit@priors$parameter,
          "a",
          fit@priors$p1,
          fit@priors$p2,
          1L,
          SimInf:::abc_particles(fit, 2L),
          fit@weight[, 2],
          sigma))
check_error(res, "Unknown distribution.")

## Check that an invalid weight is detected.
res <- assertError(
    .Call(SimInf:::SimInf_abc_proposals,
          fit@priors$parameter,
          fit@priors$distribution,
          fit@priors$p1,
          fit@priors$p2,
          1L,
          SimInf:::abc_particles(fit, 2L),
          numeric(0),
          sigma))
check_error(res, "'w' must have length >= 1 when 'x' is non-null.")

fit@weight[2, 2] <- -1
res <- assertError(
    .Call(SimInf:::SimInf_abc_proposals,
          fit@priors$parameter,
          fit@priors$distribution,
          fit@priors$p1,
          fit@priors$p2,
          1L,
          SimInf:::abc_particles(fit, 2L),
          fit@weight[, 2],
          sigma))
check_error(res, "Invalid weight detected (non-finite or < 0.0).")

fit@weight[2, 2] <- NaN
res <- assertError(
    .Call(SimInf:::SimInf_abc_proposals,
          fit@priors$parameter,
          fit@priors$distribution,
          fit@priors$p1,
          fit@priors$p2,
          1L,
          SimInf:::abc_particles(fit, 2L),
          fit@weight[, 2],
          sigma))
check_error(res, "Invalid weight detected (non-finite or < 0.0).")

fit@weight[2, 2] <- NA_real_
res <- assertError(
    .Call(SimInf:::SimInf_abc_proposals,
          fit@priors$parameter,
          fit@priors$distribution,
          fit@priors$p1,
          fit@priors$p2,
          1L,
          SimInf:::abc_particles(fit, 2L),
          fit@weight[, 2],
          sigma))
check_error(res, "Invalid weight detected (non-finite or < 0.0).")
