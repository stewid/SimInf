## siminf, a framework for stochastic disease spread simulations
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

library(siminf)

## Check measures for a SISe model
model <- SISe(init    = data.frame(id = 0, S = 99, I = 1),
              tspan   = 0:1000,
              events  = NULL,
              phi     = 1,
              upsilon = 1,
              gamma   = 0.1,
              alpha   = 1,
              beta_q1 = 1,
              beta_q2 = 1,
              beta_q3 = 1,
              beta_q4 = 1,
              epsilon = 0)

result <- run(model, threads = 1)
result

stopifnot(identical(length(susceptible(result)), 1001L))
stopifnot(identical(length(infected(result)), 1001L))
stopifnot(identical(length(prevalence(result)), 1001L))

if (siminf:::have_openmp()) {
    result_omp <- run(model, threads = 2)
    result_omp

    stopifnot(identical(length(susceptible(result_omp)), 1001L))
    stopifnot(identical(length(infected(result_omp)), 1001L))
    stopifnot(identical(length(prevalence(result_omp)), 1001L))
}

## Check measures for a SISe3 model
model <- demo_model(model = "SISe3", nodes = 10, days = 1000)

result <- run(model, threads = 1)
result

stopifnot(identical(length(susceptible(result)), 10000L))
stopifnot(identical(length(infected(result)), 10000L))

if (siminf:::have_openmp()) {
    result_omp <- run(model, threads = 2)
    result_omp

    stopifnot(identical(length(susceptible(result_omp)), 10000L))
    stopifnot(identical(length(infected(result_omp)), 10000L))
}
