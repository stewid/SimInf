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
## Create a model with parameters in ldata
##
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

accept_fn_ldata <- function(result, generation, tol, ptol, ...) {
    ## Determine the tolerance for this generation.
    tol <- tol * ptol ^ (generation - 1)

    ## Extract the time-series for R1 for each node as a
    ## data.frame.
    sim <- trajectory(result, "Icum")

    ## Split the 'sim' data.frame by node and calculate the sum of the
    ## squared distance at each time-point for every node.
    dist <- tapply(sim$Icum, sim$node, function(Icum) {
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

    ## Return TRUE or FALSE for each node depending on if the
    ## distance is less than the tolerance.
    abc_accept(dist < tol, tol)
}

## Check invalid npart
res <- assertError(abc(model = model,
                       priors = c(beta~U(0.5, 1.5), gamma~U(0.3, 0.7)),
                       ngen = 2,
                       npart = 1,
                       fn = accept_fn_ldata,
                       tol = 250000,
                       ptol = 0.9))
check_error(res, "'npart' must be an integer > 1.")

res <- assertError(abc(model = model,
                       priors = c(beta~U(0.5, 1.5), gamma~U(0.3, 0.7)),
                       ngen = 2,
                       npart = c(10, 10),
                       fn = accept_fn_ldata,
                       tol = 250000,
                       ptol = 0.9))
check_error(res, "'npart' must be an integer > 1.")

## Check invalid ngen
res <- assertError(abc(model = model,
                       priors = c(beta~U(0.5, 1.5), gamma~U(0.3, 0.7)),
                       ngen = 0,
                       npart = 10,
                       fn = accept_fn_ldata,
                       tol = 250000,
                       ptol = 0.9))
check_error(res, "'ngen' must be an integer >= 1.")

res <- assertError(abc(model = model,
                       priors = c(beta~U(0.5, 1.5), gamma~U(0.3, 0.7)),
                       ngen = c(2, 2),
                       npart = 10,
                       fn = accept_fn_ldata,
                       tol = 250000,
                       ptol = 0.9))
check_error(res, "'ngen' must be an integer >= 1.")

run_abc <- function(model) {
    ## The environmental variable "R_TEST" must be unset inside "R CMD
    ## check" in order to successfully change the working directory to
    ## a tempdir and then run "R CMD SHLIB".
    R_TESTS <- Sys.getenv("R_TESTS", unset = NA)
    if (!is.na(R_TESTS)) {
        Sys.unsetenv("R_TESTS")
        on.exit(Sys.setenv(R_TESTS = R_TESTS), add = TRUE)
    }

    abc(model = model,
        priors = c(beta~U(0.5, 1.5), gamma~U(0.3, 0.7)),
        ngen = 2,
        npart = 10,
        fn = accept_fn_ldata,
        tol = 250000,
        ptol = 0.9,
        verbose = TRUE)
}

set.seed(123)
fit <- run_abc(model)
fit
summary(fit)
as.data.frame(fit)

pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(fit, xlim = c(0.3, 1.5), ylim = c(0.3, 1.5))
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

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
