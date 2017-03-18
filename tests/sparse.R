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

U(model) <- sparseMatrix(1:18, rep(5:10, each = 3))
U_obs <- U(run(model, threads = 1, seed = 123))
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
    U(model) <- sparseMatrix(1:18, rep(5:10, each = 3))
    U_obs_omp <- U(run(model, threads = 2, seed = 123))
    stopifnot(identical(U_obs_omp, U_exp_omp))
}

## Check wrong dimension of U
m <- sparseMatrix(1:21, rep(5:11, each = 3))
res <- tools::assertError(U(model) <- m)
stopifnot(length(grep("Wrong dimension of 'value'",
                      res[[1]]$message)) > 0)

## Check wrong dimension of V
m <- as(sparseMatrix(numeric(0), numeric(0), dims = c(0, 11)), "dgCMatrix")
res <- tools::assertError(V(model) <- m)
stopifnot(length(grep("Wrong dimension of 'value'",
                      res[[1]]$message)) > 0)

## Check that U is cleared. First run a model to get a dense U result
## matrix, then run that model and check that the dense U result
## matrix is cleared. Then run the model again and check that the
## sparse result matrix is cleared.
model <- SIR(u0 = data.frame(S = 100:105, I = 1:6, R = rep(0, 6)),
             tspan = 1:10,
             beta = 0.16,
             gamma = 0.077)
result <- run(model, threads = 1)
U(result) <- sparseMatrix(1:18, rep(5:10, each = 3))
result <- run(result, threads = 1)
stopifnot(identical(dim(result@U), c(0L, 0L)))
stopifnot(identical(dim(result@U_sparse), c(18L, 10L)))
U(result) <- NULL
result <- run(result, threads = 1)
stopifnot(identical(dim(result@U), c(18L, 10L)))
stopifnot(identical(dim(result@U_sparse), c(0L, 0L)))

## Check that V is cleared. First run a model to get a dense V result
## matrix, then run that model and check that the dense V result
## matrix is cleared. Then run the model again and check that the
## sparse result matrix is cleared.
u0 <- structure(list(S  = c(0, 1, 2, 3, 4, 5),
                     I  = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S", "I"),
                row.names = c(NA, -6L), class = "data.frame")
model <- SISe(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              phi     = rep(0, nrow(u0)),
              upsilon = 0.0357,
              gamma   = 0.1,
              alpha   = 1.0,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0.000011)
result <- run(model, threads = 1)
V(result) <- sparseMatrix(1:6, 5:10)
result <- run(result, threads = 1)
stopifnot(identical(dim(result@V), c(0L, 0L)))
stopifnot(identical(dim(result@V_sparse), c(6L, 10L)))
V(result) <- NULL
result <- run(result, threads = 1)
stopifnot(identical(dim(result@V), c(6L, 10L)))
stopifnot(identical(dim(result@V_sparse), c(0L, 0L)))
