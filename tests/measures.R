## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015  Stefan Engblom
## Copyright (C) 2015  Stefan Widgren
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

## For debugging
sessionInfo()

## Check measures for a SISe model
model <- SISe(init    = data.frame(S = 99, I = 1),
              tspan   = 0:1000,
              events  = NULL,
              phi     = 1,
              upsilon = 1,
              gamma   = 0.1,
              alpha   = 1,
              beta_t1 = 1,
              beta_t2 = 1,
              beta_t3 = 1,
              beta_t4 = 1,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0)

res <- tools::assertError(susceptible(model))
stopifnot(length(grep("Please run the model first, the 'U' matrix is empty",
                      res[[1]]$message)) > 0)

res <- tools::assertError(infected(model))
stopifnot(length(grep("Please run the model first, the 'U' matrix is empty",
                      res[[1]]$message)) > 0)

result <- run(model, threads = 1)
result

stopifnot(identical(length(susceptible(result)), 1001L))
i <- seq(from = 1, to = 1001, by = 2)
stopifnot(identical(
    susceptible(result)[, i, drop = FALSE],
    susceptible(result, by = 2)))

stopifnot(identical(length(infected(result)), 1001L))
stopifnot(identical(
    infected(result)[, i, drop = FALSE],
    infected(result, by = 2)))

stopifnot(identical(length(prevalence(result)), 1001L))
stopifnot(is.null(dim(prevalence(result))))
stopifnot(identical(dim(prevalence(result, wnp = TRUE)), c(1L, 1001L)))
stopifnot(identical(
    prevalence(result)[i],
    prevalence(result, by = 2)))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2)
    result_omp

    stopifnot(identical(length(susceptible(result_omp)), 1001L))
    i <- seq(from = 1, to = 1001, by = 2)
    stopifnot(identical(
        susceptible(result_omp)[, i, drop = FALSE],
        susceptible(result_omp, by = 2)))

    stopifnot(identical(length(infected(result_omp)), 1001L))
    stopifnot(identical(
        infected(result_omp)[, i, drop = FALSE],
        infected(result_omp, by = 2)))

    stopifnot(identical(length(prevalence(result_omp)), 1001L))
    stopifnot(is.null(dim(prevalence(result_omp))))
    stopifnot(identical(dim(prevalence(result_omp, wnp = TRUE)), c(1L, 1001L)))
    stopifnot(identical(
        prevalence(result_omp)[i],
        prevalence(result_omp, by = 2)))
}

## Check measures for a SISe3 model
model <- demo_model(model = "SISe3", nodes = 10, days = 1000)

res <- tools::assertError(susceptible(model))
stopifnot(length(grep("Please run the model first, the 'U' matrix is empty",
                      res[[1]]$message)) > 0)

res <- tools::assertError(infected(model))
stopifnot(length(grep("Please run the model first, the 'U' matrix is empty",
                      res[[1]]$message)) > 0)

result <- run(model, threads = 1)
result

stopifnot(identical(length(susceptible(result)), 10000L))
i <- seq(from = 1, to = 1000, by = 2)
stopifnot(identical(
    susceptible(result, age = 1)[, i, drop = FALSE],
    susceptible(result, age = 1, by = 2)))

stopifnot(identical(length(infected(result)), 10000L))
stopifnot(identical(
    infected(result, age = 1)[, i, drop = FALSE],
    infected(result, age = 1, by = 2)))

stopifnot(identical(length(prevalence(result)), 1000L))
stopifnot(is.null(dim(prevalence(result))))
stopifnot(identical(dim(prevalence(result, wnp = TRUE)), c(10L, 1000L)))
stopifnot(identical(
    prevalence(result)[i],
    prevalence(result, by = 2)))

if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2)
    result_omp

    stopifnot(identical(length(susceptible(result_omp)), 10000L))
    i <- seq(from = 1, to = 1000, by = 2)
    stopifnot(identical(
        susceptible(result_omp)[, i, drop = FALSE],
        susceptible(result_omp, by = 2)))

    stopifnot(identical(length(infected(result_omp)), 10000L))
    stopifnot(identical(
        infected(result_omp)[, i, drop = FALSE],
        infected(result_omp, by = 2)))

    stopifnot(identical(length(prevalence(result_omp)), 1000L))
    stopifnot(is.null(dim(prevalence(result_omp))))
    stopifnot(identical(dim(prevalence(result_omp, wnp = TRUE)), c(10L, 1000L)))
    stopifnot(identical(
        prevalence(result_omp)[i],
        prevalence(result_omp, by = 2)))
}
