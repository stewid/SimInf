## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2017  Robin Eriksson
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

## For debugging
sessionInfo()

## Check measures with a SIR model
u0 <- structure(list(S  = c(0, 1, 2, 3, 4, 5),
                     I  = c(0, 0, 0, 0, 0, 0),
                     R  = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S", "I", "R"),
                row.names = c(NA, -6L), class = "data.frame")

model <- SIR(u0     = u0,
             tspan  = seq_len(10) - 1,
             events = NULL,
             beta   = 0.16,
             gamma  = 0.077)

res <- tools::assertError(susceptible(model))
stopifnot(length(grep("Please run the model first, the 'U' matrix is empty",
                      res[[1]]$message)) > 0)

res <- tools::assertError(infected(model))
stopifnot(length(grep("Please run the model first, the 'U' matrix is empty",
                      res[[1]]$message)) > 0)

## run with AEM
result <- run(model, threads = 1, solver = "aem")
result

stopifnot(identical(length(susceptible(result)), 60L))
stopifnot(identical(length(infected(result)), 60L))
stopifnot(identical(length(prevalence(result)), 10L))
stopifnot(is.null(dim(prevalence(result))))
stopifnot(identical(dim(prevalence(result, type = "wnp")), c(6L, 10L)))

## run with AEM using multiple threads
if (SimInf:::have_openmp()) {
    result_omp <- run(model, threads = 2, solver = "aem")
    result_omp

    stopifnot(identical(length(susceptible(result_omp)), 60L))
    stopifnot(identical(length(infected(result_omp)), 60L))
    stopifnot(identical(length(prevalence(result_omp)), 10L))
    stopifnot(is.null(dim(prevalence(result_omp))))
    stopifnot(identical(dim(prevalence(result_omp, type = "wnp")), c(6L, 10L)))
}
