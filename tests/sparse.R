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

library("SimInf")
library("Matrix")

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
                             c("1", "2", "3", "4", "5",
                               "6", "7", "8", "9", "10")),
             x = c(100, 0, 1, 101, 2, 0, 102, 0, 3, 98, 7, 2,
                   88, 11, 10, 101, 5, 5),
             factors = list())

U(model) <- sparseMatrix(1:18, rep(5:10, each = 3))
U_obs <- trajectory(run(model, threads = 1, seed = 123), as.is = TRUE)
stopifnot(identical(U_obs, U_exp))

if (SimInf:::have_openmp()) {
    U_exp_omp <- new("dgCMatrix",
                     i = 0:17,
                     p = c(0L, 0L, 0L, 0L, 0L, 3L, 6L, 9L, 12L, 15L, 18L),
                     Dim = c(18L, 10L),
                     Dimnames = list(c("S", "I", "R", "S", "I", "R",
                                       "S", "I", "R", "S", "I", "R",
                                       "S", "I", "R", "S", "I", "R"),
                                     c("1", "2", "3", "4", "5",
                                       "6", "7", "8", "9", "10")),
                     x = c(96, 5, 0, 101, 1, 1, 102, 1, 2,
                           99, 6, 2, 98, 3, 8, 95, 12, 4),
                     factors = list())
    U(model) <- sparseMatrix(1:18, rep(5:10, each = 3))
    U_obs_omp <- trajectory(run(model, threads = 2, seed = 123), as.is = TRUE)
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

## Check data.frame output from sparse matrix. Create an 'SIR' model
## with 6 nodes and initialize it to run over 10 days. Then create a
## sparse matrix with non-zero entries at the locations in U where the
## number of individuals should be written. Run the model with the
## sparse matrix as a template for U where to write data.
u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
m <- Matrix::sparseMatrix(1:12, rep(5:10, each = 2), dims = c(18, 10))
U(model) <- m
result <- run(model, threads = 1, seed = 22)
U_exp <- structure(list(Node = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L),
                        Time = c(5L, 6L, 6L, 7L, 8L, 9L, 9L, 10L),
                        S = c(98L, NA, 100L, NA, 96L, NA, 101L, NA),
                        I = c(3L, NA, NA, 3L, 7L, NA, NA, 3L),
                        R = c(NA, 0L, NA, 1L, NA, 3L, NA, 3L)),
                   .Names = c("Node", "Time", "S", "I", "R"),
                   row.names = c(NA, -8L),
                   class = "data.frame")
stopifnot(identical(trajectory(result), U_exp))

## Similar test case, but without NA-values
m <- Matrix::sparseMatrix(1:18, rep(5:10, each = 3))
U(model) <- m
result <- run(model, threads = 1, seed = 22)
U_exp <- structure(list(Node = 1:6,
                        Time = 5:10,
                        S = c(98L, 100L, 97L, 101L, 87L, 93L),
                        I = c(3L, 2L, 6L, 3L, 15L, 16L),
                        R = c(0L, 1L, 2L, 3L, 7L, 2L)),
                   .Names = c("Node", "Time", "S", "I", "R"),
                   row.names = c(NA, -6L),
                   class = "data.frame")
stopifnot(identical(trajectory(result), U_exp))

## Test that it also works to remove the sparse matrix output
U(model) <- NULL
result <- run(model, threads = 1, seed = 22)
U_exp <- structure(list(
    Node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L,
             3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L,
             1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,  5L, 6L),
    Time = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L,
             3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L,
             6L, 6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 8L,
             9L, 9L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 10L, 10L, 10L),
    S = c(100L, 101L, 102L, 103L, 104L, 105L, 99L, 101L, 102L, 102L, 102L,
          105L, 98L, 101L, 102L, 101L, 100L, 105L, 98L, 101L, 99L, 101L, 97L,
          104L, 98L, 100L, 99L, 101L, 97L, 102L, 98L, 100L, 98L, 101L, 95L,
          99L, 98L, 99L, 97L, 101L, 93L, 98L, 98L, 97L, 96L, 101L, 92L, 95L,
          97L, 96L, 95L, 101L, 87L, 94L, 97L, 95L, 95L, 101L, 86L, 93L),
    I = c(1L, 2L, 3L, 4L, 5L, 6L, 2L, 2L, 3L, 5L, 7L, 6L, 3L, 2L, 3L, 6L, 8L,
          6L, 3L, 1L, 4L, 3L, 8L, 7L, 3L, 2L, 4L, 3L, 8L, 9L, 3L, 2L, 5L, 3L,
          8L, 11L, 3L, 3L, 6L, 3L, 9L, 12L, 2L, 5L, 7L, 3L, 10L, 15L, 3L, 6L,
          7L, 3L, 15L, 15L, 3L, 6L, 4L, 3L, 13L, 16L),
    R = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L,
          0L, 0L, 1L, 2L, 3L, 4L, 0L, 0L, 1L, 2L, 3L, 4L, 0L, 0L, 1L, 2L, 3L,
          6L, 1L, 0L, 1L, 2L, 3L, 7L, 1L, 1L, 1L, 2L, 3L, 7L, 1L, 1L, 1L, 3L,
          3L, 7L, 2L, 1L, 2L, 6L, 3L, 10L, 2L)),
    .Names = c("Node", "Time", "S", "I", "R"),
    row.names = c(NA, -60L),
    class = "data.frame")
stopifnot(identical(trajectory(result), U_exp))
