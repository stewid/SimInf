## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2021 Stefan Widgren
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

res <- assertError(SimInf:::parse_priors(c(a + b ~ uniform(0, 5),
                                           c ~ uniform(0, 1), 4)))
check_error(res, "'priors' must be a formula or a list with formula items.")

res <- assertError(SimInf:::parse_priors(4))
check_error(res, "'priors' must be a formula or a list with formula items.")

res <- assertError(SimInf:::parse_priors(NULL))
check_error(res, "'priors' must be a formula or a list with formula items.")

res <- assertError(SimInf:::parse_priors(mu ~ uniform(0, 1) + normal(0, 1)))
check_error(res, "Invalid formula specification for distribution.")

res <- assertError(SimInf:::parse_priors(mu ~ uniform[0, 1]))
check_error(res, "Invalid formula specification for distribution.")

res <- assertError(SimInf:::parse_priors(mu ~ unknown(0, 1)))
check_error(
    res, "'distribution' must be one of 'gamma', 'normal' or 'uniform'.")

res <- assertError(SimInf:::parse_priors(c(muR ~ uniform(0, 1),
                                           muR ~ uniform(0, 1))))
check_error(res, "'priors' must have non-duplicated parameter names.")

res <- assertError(SimInf:::parse_priors(beta ~ uniform(1, 0)))
check_error(res, "Invalid distribution: uniform bounds in wrong order.")

res <- assertError(SimInf:::parse_priors(beta ~ normal(0, -1)))
check_error(res, "Invalid distribution: normal variance must be > 0.")

res <- assertError(SimInf:::parse_priors(beta ~ gamma(-1, 1)))
check_error(res, "Invalid distribution: gamma hyperparameters must be > 0.")

res <- assertError(SimInf:::parse_priors(beta ~ gamma(1, -1)))
check_error(res, "Invalid distribution: gamma hyperparameters must be > 0.")

res <- assertError(SimInf:::parse_priors(~ uniform(1, 5)))
check_error(res, "Invalid formula specification for distribution.")

stopifnot(identical(
    SimInf:::parse_priors(beta ~ uniform(1, 5)),
    data.frame(parameter = "beta", distribution = "uniform",
               p1 = 1, p2 = 5, stringsAsFactors = FALSE)))
