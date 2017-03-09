## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2017  Stefan Engblom
## Copyright (C) 2015 - 2017  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(SimInf)
library(Matrix)

## For debugging
sessionInfo()

## Check sparse U
model <- SIR(u0 = data.frame(S = 100:105, I = 1:6, R = rep(0, 6)),
             tspan = 1:10,
             beta = 0.16,
             gamma = 0.077)

U_exp <- new("dgCMatrix",
             i = 0:17,
             p = c(0L, 0L, 0L, 0L, 0L, 3L, 6L, 9L, 12L, 15L, 18L),
             Dim = c(18L, 10L),
             Dimnames = list(c("S", "I", "R", "S", "I", "R",
                               "S", "I", "R", "S", "I", "R",
                               "S", "I", "R", "S", "I", "R"),
                             NULL),
             x = c(100, 0, 1, 100, 3, 0, 97, 7, 1,
                   94, 10, 3, 91, 13, 5, 92, 12, 7),
             factors = list())

m <- sparseMatrix(1:18, rep(5:10, each = 3))
U_obs <- U(run(model, threads = 1, seed = 123, U = m))
stopifnot(identical(U_obs, U_exp))

if (SimInf:::have_openmp()) {
    U_exp_omp <- new("dgCMatrix",
                     i = 0:17,
                     p = c(0L, 0L, 0L, 0L, 0L, 3L, 6L, 9L, 12L, 15L, 18L),
                     Dim = c(18L, 10L),
                     Dimnames = list(c("S", "I", "R", "S", "I", "R",
                                       "S", "I", "R", "S", "I", "R",
                                       "S", "I", "R", "S", "I", "R"),
                                     NULL),
                     x = c(97, 4, 0, 100, 3, 0, 97, 7, 1,
                           92, 9, 6, 98, 7, 4, 94, 9, 8),
                     factors = list())
    U_obs_omp <- U(run(model, threads = 2, seed = 123, U = m))
    stopifnot(identical(U_obs_omp, U_exp_omp))
}

## Check wrong dimension of U
m <- sparseMatrix(1:21, rep(5:11, each = 3))
res <- tools::assertError(run(model, U = m))
stopifnot(length(grep("Wrong dimension of 'U'",
                      res[[1]]$message)) > 0)

## Check wrong dimension of V
m <- as(sparseMatrix(numeric(0), numeric(0), dims = c(0, 11)), "dgCMatrix")
res <- tools::assertError(run(model, V = m))
stopifnot(length(grep("Wrong dimension of 'V'",
                      res[[1]]$message)) > 0)
