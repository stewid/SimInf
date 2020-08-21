## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
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

library(SimInf)
library(tools)
source("util/check.R")

## Specify the number of threads to use.
set_num_threads(1)

##
## Create a model with parameters in gdata
##
model <- mparse(transitions = c("S -> beta*S*I/(S+I+R) -> I",
                                "I -> gamma*I -> R"),
                compartments = c("S", "I", "R"),
                gdata = c(beta = 1, gamma = 0.5),
                u0 = data.frame(S = rep(9999, 2), I = 1, R = 0),
                tspan = 1:50)

accept_fn_gdata <- function(result, generation, tol, ptol, ...) {
    ## Determine the tolerance for this generation.
    tol <- tol * ptol ^ (generation - 1)

    p <- c(2e-04, 0.00015, 5e-05, 5e-05, 2e-04, 0.00025, 0.00025,
           0.00025, 0.00025, 0.00015, 0.00035, 6e-04, 0.001, 0.0022,
           0.00395, 0.00655, 0.0102, 0.01755, 0.02795, 0.04235,
           0.05925, 0.07135, 0.08025, 0.08205, 0.0744, 0.0657,
           0.05785, 0.04775, 0.03735, 0.02855, 0.02265, 0.01775,
           0.0128, 0.01005, 0.00745, 0.00545, 0.0038, 0.0027, 0.00205,
           0.00145, 0.0012, 8e-04, 7e-04, 3e-04, 2e-04, 0.00015,
           5e-05, 0, 0, 0)

    dist <- sum((prevalence(result, I~.)$prevalence - p)^2)

    ## Return TRUE or FALSE depending on if the distance is less than
    ## or equal to the tolerance.
    abc_accept(dist < tol, tol)
}

set.seed(123)
fit <- abc(model = model,
           priors = c(beta~U(0.5, 1.5), gamma~U(0.3, 0.7)),
           ngen = 2,
           npart = 10,
           fn = accept_fn_gdata,
           tol = 0.1,
           ptol = 0.5)
fit
summary(fit)

## Test that PKG_CPPFLAGS is unset.
stopifnot(isTRUE(is.na(Sys.getenv("PKG_CPPFLAGS", unset = NA))))
