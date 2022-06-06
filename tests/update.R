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

model <- SIR(u0     = data.frame(S = 1:3, I = 4:6, R = 7:9),
             tspan  = 1:10,
             events = NULL,
             beta   = 0,
             gamma  = 0)

model <- update_u0(model, data.frame(S = 10:12, I = 13:15, R = 16:18))
stopifnot(identical(model@u0,
                    matrix(c(10L, 11L, 12L,
                             13L, 14L, 15L,
                             16L, 17L, 18L),
                           nrow = 3,
                           ncol = 3,
                           byrow = TRUE,
                           dimnames = list(c("S", "I", "R"), NULL))))

res <- assertError(
    update_u0(model, data.frame(S = 10:13, I = 14:17, R = 18:21)))
check_error(res, "The number of rows in 'u0' must match nodes in 'model'.")

res <- assertError(
    update_v0(model, data.frame(phi = 10:13)))
check_error(res, "The number of rows in 'v0' must match nodes in 'model'.")

model <- update_v0(model, data.frame(phi = 1:3))
stopifnot(identical(model@v0, matrix(numeric(0), nrow = 0, ncol = 0)))

## Check that a matrix is coerced to a data.frame.
res <- SimInf:::check_v0(matrix(1:9,
                                ncol = 3,
                                dimnames = list(NULL, c("A", "B", "C"))),
                         c("A", "B", "C"))
stopifnot(identical(res,
                    data.frame(A = 1:3,
                               B = 4:6,
                               C = 7:9)))

## Create an 'SISe' model with 6 nodes.
model <- SISe(u0 = data.frame(S = 100:105, I = 1:6), tspan = 1:10,
              phi = rep(0, 6), upsilon = 0.02, gamma = 0.1, alpha = 1,
              epsilon = 1.1e-5, beta_t1 = 0.15, beta_t2 = 0.15,
              beta_t3 = 0.15, beta_t4 = 0.15, end_t1 = 91, end_t2 = 182,
              end_t3 = 273, end_t4 = 365)

res <- assertError(update_v0(model, data.frame(A = 1:6)))
check_error(res, "Missing columns in 'v0'.")

model <- update_v0(model, data.frame(phi = 1:6))
stopifnot(identical(model@v0,
                    matrix(c(1, 2, 3, 4, 5, 6),
                           nrow = 1,
                           ncol = 6,
                           dimnames = list("phi", NULL))))
