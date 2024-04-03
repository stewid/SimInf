## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2024 Stefan Widgren
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

## For debugging
sessionInfo()

model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = seq(1, 100, by = 3),
             beta = 0.16,
             gamma = 0.077)

## Observed data
infected <- data.frame(
    time = c(1L, 4L, 7L, 10L, 13L, 16L, 19L, 22L, 25L, 28L, 31L, 34L,
             37L, 40L, 43L, 46L, 49L, 52L, 55L, 58L, 61L, 64L, 67L,
             70L, 73L, 76L, 79L, 82L, 85L, 88L, 91L, 94L, 97L, 100L),
    Iobs = c(1L, 2L, 2L, 3L, 3L, 3L, 3L, 2L, 3L, 3L, 6L, 10L, 9L, 16L,
             17L, 15L, 15L, 13L, 18L, 18L, 16L, 13L, 12L, 12L, 11L,
             10L, 9L, 10L, 7L, 7L, 6L, 3L, 3L, 2L))

## Check that an invalid 'npart' raises an error.
res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          npart = 200.1,
          niter = 200,
          verbose = TRUE))
check_error(res, "'npart' must be integer.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          npart = NA_integer_,
          niter = 200,
          verbose = TRUE))
check_error(res, "'npart' must be integer.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          npart = -200,
          niter = 200,
          verbose = TRUE))
check_error(res, "'npart' must be an integer > 1.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          npart = c(200, 200),
          niter = 200,
          verbose = TRUE))
check_error(res, "'npart' must be an integer > 1.")

## Check that an invalid 'adaptmix' raises an error.
res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          npart = 200,
          niter = 200,
          adaptmix = -1,
          verbose = TRUE))
check_error(res, "'adaptmix' must be a value > 0 and < 1.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          npart = 200,
          niter = 200,
          adaptmix = c(0.5, 0.5),
          verbose = TRUE))
check_error(res, "'adaptmix' must be a value > 0 and < 1.")

## Check that an invalid 'niter' raises an error.
res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          npart = 200,
          niter = 0,
          theta = c(beta = 0.16, gamma = 0.077)))
check_error(res, "'niter' must be an integer > 0.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          npart = 200,
          niter = c(200, 200),
          theta = c(beta = 0.16, gamma = 0.077)))
check_error(res, "'niter' must be an integer > 0.")

## Check that an invalid 'theta' raises an error.
res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          npart = 200,
          niter = 200,
          theta = c(beta = "A", gamma = 0.077)))
check_error(
    res,
    "'theta' must be a vector with initial values for the parameters.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          npart = 200,
          niter = 200,
          theta = c(gamma = 0.077)))
check_error(
    res,
    "'theta' must be a vector with initial values for the parameters.")

## Run pmcmc
set.seed(123)
fit <- pmcmc(model,
             Iobs ~ poisson(I + 1e-6),
             infected,
             priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
             npart = 10,
             niter = 1,
             theta = c(beta = 0.16, gamma = 0.077))

show_expected <- c(
    "Particle Markov chain Monte Carlo",
    "---------------------------------",
    "Number of iterations: 1",
    "Number of particles: 10",
    "Mixing proportion for adaptive proposal: 0.05",
    "Acceptance ratio: 0.000",
    "",
    "Quantiles, mean and standard deviation for each variable",
    "--------------------------------------------------------",
    "       2.5%   25%   50%   75% 97.5%  Mean SD",
    "beta  0.160 0.160 0.160 0.160 0.160 0.160   ",
    "gamma 0.077 0.077 0.077 0.077 0.077 0.077   ")
show_observed <- capture.output(show(fit))
stopifnot(identical(show_observed, show_expected))

stopifnot(isTRUE(SimInf:::valid_SimInf_pmcmc_object(fit)))

fit@adaptmix <- 1:2
stopifnot(identical(SimInf:::valid_SimInf_pmcmc_object(fit),
                    "'adaptmix' must be a value >= 0 and <= 1."))

fit@adaptmix <- 0
stopifnot(identical(SimInf:::valid_SimInf_pmcmc_object(fit),
                    "'adaptmix' must be a value >= 0 and <= 1."))

fit@adaptmix <- 1
stopifnot(identical(SimInf:::valid_SimInf_pmcmc_object(fit),
                    "'adaptmix' must be a value >= 0 and <= 1."))

fit@adaptmix <- 0.05
fit@target <- "test"
stopifnot(identical(SimInf:::valid_SimInf_pmcmc_object(fit),
                    "'target' must be 'gdata' or 'ldata'."))

fit <- pmcmc(model,
             Iobs ~ poisson(I + 1e-6),
             infected,
             priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
             npart = 10,
             niter = 5)
summary_expected <- c(
    "Particle Markov chain Monte Carlo",
    "---------------------------------",
    "Number of iterations: 5",
    "Number of particles: 10",
    "Acceptance ratio: 0.400",
    "Model: SIR", "Number of nodes: 1",
    "",
    "Transitions",
    "-----------",
    " S -> beta*S*I/(S+I+R) -> I",
    " I -> gamma*I -> R",
    "",
    "Quantiles, mean and standard deviation for each variable",
    "--------------------------------------------------------",
    "        2.5%    25%    50%    75%  97.5%   Mean     SD",
    "beta  0.7545 0.7545 0.8035 0.8663 0.8663 0.8090 0.0560",
    "gamma 0.6062 0.6292 0.6292 0.6438 0.6438 0.6299 0.0164")
summary_observed <- capture.output(summary(fit))
stopifnot(identical(summary_observed, summary_expected))

res <- assertError(
    continue(fit, niter = 0))
check_error(
    res,
    "'niter' must be an integer > 0.")

res <- assertError(
    continue(fit, niter = 1:2))
check_error(
    res,
    "'niter' must be an integer > 0.")

stopifnot(identical(
    SimInf:::pmcmc_iterations(x = fit, start = 1, end = NULL, thin = 1),
    1:5))

res <- assertError(
    SimInf:::pmcmc_iterations(x = fit, start = -1, end = NULL, thin = 1))
check_error(
    res,
    "'start' must be an integer >= 1.")

res <- assertError(
    SimInf:::pmcmc_iterations(x = fit, start = 1:2, end = NULL, thin = 1))
check_error(
    res,
    "'start' must be an integer >= 1.")

res <- assertError(
    SimInf:::pmcmc_iterations(x = fit, start = 1, end = 1:2, thin = 1))
check_error(
    res,
    "'end' must be an integer between start and length(x).")

res <- assertError(
    SimInf:::pmcmc_iterations(x = fit, start = 2, end = 1, thin = 1))
check_error(
    res,
    "'end' must be an integer between start and length(x).")

res <- assertError(
    SimInf:::pmcmc_iterations(x = fit, start = 1, end = 6, thin = 1))
check_error(
    res,
    "'end' must be an integer between start and length(x).")

res <- assertError(
    SimInf:::pmcmc_iterations(x = fit, start = 1, end = NULL, thin = 1:2))
check_error(
    res,
    "'thin' must be an integer >= 1.")

res <- assertError(
    SimInf:::pmcmc_iterations(x = fit, start = 1, end = NULL, thin = -1))
check_error(
    res,
    "'thin' must be an integer >= 1.")
