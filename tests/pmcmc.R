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

## Define a tolerance
tol <- 1e-8

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

## Check that an invalid 'n_particles' raises an error.
res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          n_particles = 200.1,
          niter = 200,
          verbose = TRUE))
check_error(res, "'n_particles' must be integer.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          n_particles = NA_integer_,
          niter = 200,
          verbose = TRUE))
check_error(res, "'n_particles' must be integer.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          n_particles = -200,
          niter = 200,
          verbose = TRUE))
check_error(res, "'n_particles' must be an integer > 1.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          n_particles = c(200, 200),
          niter = 200,
          verbose = TRUE))
check_error(res, "'n_particles' must be an integer > 1.")

## Check that an invalid 'adaptmix' raises an error.
res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          n_particles = 200,
          niter = 200,
          adaptmix = -1,
          verbose = TRUE))
check_error(res, "'adaptmix' must be a value > 0 and < 1.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          n_particles = 200,
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
          n_particles = 200,
          niter = 0,
          theta = c(beta = 0.16, gamma = 0.077)))
check_error(res, "'n_iterations' must be an integer > 0.")

res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          n_particles = 200,
          niter = c(200, 200),
          theta = c(beta = 0.16, gamma = 0.077)))
check_error(res, "'n_iterations' must be an integer > 0.")

## Check that an invalid 'theta' raises an error.
res <- assertError(
    pmcmc(model,
          Iobs ~ poisson(I + 1e-6),
          infected,
          priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
          n_particles = 200,
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
          n_particles = 200,
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
             n_particles = 10,
             niter = 5,
             theta = c(beta = 0.16, gamma = 0.077))

show_expected <- c(
    "Particle Markov chain Monte Carlo",
    "---------------------------------",
    "Number of iterations: 5",
    "Number of particles: 10",
    "Mixing proportion for adaptive proposal: 0.05",
    "Acceptance ratio: 0.400",
    "",
    "Quantiles, mean and standard deviation for each variable",
    "--------------------------------------------------------",
    "         2.5%     25%     50%     75%   97.5%    Mean      SD",
    "beta  0.15252 0.16000 0.16000 0.16749 0.16749 0.16133 0.00657",
    "gamma 0.07000 0.07700 0.07700 0.07883 0.07883 0.07618 0.00399")
show_observed <- capture.output(show(fit))
stopifnot(identical(show_observed, show_expected))

summary_expected <- c(
    "Particle Markov chain Monte Carlo",
    "---------------------------------",
    "Number of iterations: 5",
    "Number of particles: 10",
    "Acceptance ratio: 0.400",
    "Model: SIR",
    "Number of nodes: 1",
    "",
    "Transitions",
    "-----------",
    " S -> beta*S*I/(S+I+R) -> I",
    " I -> gamma*I -> R",
    "",
    "Quantiles, mean and standard deviation for each variable",
    "--------------------------------------------------------",
    "         2.5%     25%     50%     75%   97.5%    Mean      SD",
    "beta  0.15252 0.16000 0.16000 0.16749 0.16749 0.16133 0.00657",
    "gamma 0.07000 0.07700 0.07700 0.07883 0.07883 0.07618 0.00399")
summary_observed <- capture.output(summary(fit))
stopifnot(identical(summary_observed, summary_expected))

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

stopifnot(all(is.na(SimInf:::setup_chain(fit, 5)[6:10, ])))

fit@target <- "ldata"
res <- assertError(
    continue_pmcmc(fit, niter = 0))
check_error(
    res,
    "'n_iterations' must be an integer > 0.")

res <- assertError(
    continue_pmcmc(fit, niter = 1:2))
check_error(
    res,
    "'n_iterations' must be an integer > 0.")

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

set.seed(22)
fit@chain <- fit@chain[sample(1:5, 100, replace = TRUE), ]
progress_expected <- c(
"             2.5%       25%       50%       75%     97.5%      Mean        SD",
"logPost -79.76966 -79.76966 -79.45949 -79.45949 -78.50691 -79.37731   0.47179",
"beta      0.15168   0.16000   0.16000   0.16749   0.16749   0.16110   0.00592",
"gamma     0.06923   0.07700   0.07700   0.07883   0.07883   0.07606   0.00364")
progress_observed <- capture.output(SimInf:::pmcmc_progress(fit, 100, TRUE))
## Skip first three lines since it contains a timestamp.
progress_observed <- progress_observed[4:7]
stopifnot(identical(progress_observed, progress_expected))

plot(fit)
plot(fit, ~trace)

fit@model@gdata <- c(beta = 0, gamma = 0)
fit@target <- "gdata"
stopifnot(all(
    abs(SimInf:::set_proposal(fit, c(beta = 0.5, gamma = 0.6)) -
        c(beta = 0.5, gamma = 0.6)) < tol))
