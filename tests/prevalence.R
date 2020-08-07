## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
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

## For debugging
sessionInfo()

## Define a tolerance
tol <- 1e-8

model <- SIR(u0 = data.frame(S = c(8, 5, 0), I = c(0, 1, 0), R = c(0, 0, 4)),
             tspan = 1:5, beta = 0.1, gamma = 0.1)

res <- assertError(prevalence(model, I ~ . | R == 0))
check_error(res, "Please run the model first, the trajectory is empty.")

model@U <- matrix(c(8L, 8L, 8L, 8L, 8L,
                    0L, 0L, 0L, 0L, 0L,
                    0L, 0L, 0L, 0L, 0L,
                    5L, 4L, 3L, 2L, 1L,
                    1L, 2L, 3L, 3L, 3L,
                    0L, 0L, 0L, 1L, 2L,
                    0L, 0L, 0L, 0L, 0L,
                    0L, 0L, 0L, 0L, 0L,
                    4L, 4L, 4L, 4L, 4L),
                  ncol = 5,
                  byrow = TRUE,
                  dimnames = list(c("S", "I", "R",
                                    "S", "I", "R",
                                    "S", "I", "R"),
                                  c("1", "2", "3", "4", "5")))

res <- assertError(prevalence(model, ~I))
check_error(res, "Invalid 'formula' specification.")

p <- prevalence(model, I ~ .)$prevalence
stopifnot(all(abs(p - c(1 / 18, 2 / 18, 3 / 18, 3 / 18, 3 / 18)) < tol))

p <- prevalence(model, I ~ . | R == 0)$prevalence
stopifnot(all(abs(p - c(1 / 14, 2 / 14, 3 / 14, 0 / 8, 0 / 8)) < tol))

p <- prevalence(model, I ~ . | R > 0)$prevalence
stopifnot(all(abs(p - c(0 / 4, 0 / 4, 0 / 4, 3 / 10, 3 / 10)) < tol))

stopifnot(all(is.nan(prevalence(model, I ~ . | R == 5)$prevalence)))

res <- assertError(prevalence(model, I ~ . | TRUE == 0))
check_error(
    res,
    paste("The condition must be either 'TRUE' or",
          "'FALSE' for every node and time step."))

p <- prevalence(model, I ~ . | S == 0 | R == 0)$prevalence
stopifnot(all(abs(p - c(1 / 18, 2 / 18, 3 / 18, 0 / 12, 0 / 12)) < tol))

p <- prevalence(model, I ~ . | S == 0 | R == 0, i = 2)$prevalence
stopifnot(all(abs(p[1:3] - c(1 / 6, 2 / 6, 3 / 6)) < tol))
stopifnot(all(is.nan(p[4:5])))
