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

stopifnot(all(is.na(SimInf:::setup_chain(fit, 5)[6:10, ])))

proposal_exp <- list(
    theta = c(beta = 0.805085162476816, gamma = 0.672626293669262),
    theta_mean = c(beta = 0.818561741932892, gamma = 0.632232762500019),
    covmat_emp = structure(c(0.00299324946584347, 0.000468732158843136,
                             0.000468732158843136, 0.000246502076160179),
                           dim = c(2L, 2L),
                           dimnames = list(c("beta", "gamma"),
                                           c("beta", "gamma"))))
theta_mean <- colMeans(fit@chain[seq_len(5), 5:6, drop = FALSE])
covmat_emp <- SimInf:::covmat_empirical(fit, 5)
proposal_obs <- SimInf:::pmcmc_proposal(fit, i = 6, n_accepted = 2,
                                        theta_mean = theta_mean,
                                        covmat_emp = covmat_emp,
                                        scale_start = 5,
                                        shape_start = 200,
                                        scale_cooling = 0.999,
                                        max_scaling = 50)

stopifnot(all(abs(proposal_exp$theta - proposal_obs$theta) < tol))
stopifnot(all(abs(proposal_exp$theta_mean - proposal_obs$theta_mean) < tol))
stopifnot(all(abs(proposal_exp$covmat_emp - proposal_obs$covmat_emp) < tol))

proposal_exp <- list(
    theta = c(beta = 0.832033430879964, gamma = 0.640270855570007),
    theta_mean = c(beta = 0.818561741932892, gamma = 0.632232762500019),
    covmat_emp = structure(c(0.00299324946584347, 0.000468732158843136,
                             0.000468732158843136, 0.000246502076160179),
                           dim = c(2L, 2L),
                           dimnames = list(c("beta", "gamma"),
                                           c("beta", "gamma"))))
theta_mean <- colMeans(fit@chain[seq_len(5), 5:6, drop = FALSE])
covmat_emp <- SimInf:::covmat_empirical(fit, 5)
proposal_obs <- SimInf:::pmcmc_proposal(fit, i = 6, n_accepted = 2,
                                        theta_mean = theta_mean,
                                        covmat_emp = covmat_emp,
                                        scale_start = 5,
                                        shape_start = 2,
                                        scale_cooling = 0.999,
                                        max_scaling = 50)

stopifnot(all(abs(proposal_exp$theta - proposal_obs$theta) < tol))
stopifnot(all(abs(proposal_exp$theta_mean - proposal_obs$theta_mean) < tol))
stopifnot(all(abs(proposal_exp$covmat_emp - proposal_obs$covmat_emp) < tol))

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

set.seed(22)
fit@chain <- fit@chain[sample(1:5, 100, replace = TRUE), ]
progress_expected <- c(
"",
"PMCMC iteration: 100 of 100. Acceptance ratio: 0.360",
"----------------------------------------------------",
"             2.5%       25%       50%       75%     97.5%      Mean        SD",
"logPost -3.94e+03 -3.94e+03 -3.68e+03 -3.20e+03 -3.20e+03 -3.60e+03  3.36e+02",
"beta     7.54e-01  7.54e-01  8.04e-01  8.66e-01  8.66e-01  8.07e-01  5.00e-02",
"gamma    6.04e-01  6.29e-01  6.29e-01  6.44e-01  6.44e-01  6.29e-01  1.48e-02")
progress_observed <- capture.output(SimInf:::pmcmc_progress(fit, 100, TRUE))
stopifnot(identical(progress_observed, progress_expected))

plot(fit)
plot(fit, ~trace)

fit@model@gdata <- c(beta = 0, gamma = 0)
fit@target <- "gdata"
stopifnot(all(
    abs(SimInf:::set_proposal(fit, c(beta = 0.5, gamma = 0.6)) -
        c(beta = 0.5, gamma = 0.6)) < tol))
