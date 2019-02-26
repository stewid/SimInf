## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015 - 2019  Stefan Engblom
## Copyright (C) 2015 - 2019  Stefan Widgren
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
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

library("SimInf")

## For debugging
sessionInfo()

## Check measures for a SISe model
model <- SISe(u0      = data.frame(S = 99, I = 1),
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

res <- tools::assertError(trajectory(model, compartments = "S", as.is = TRUE))
stopifnot(length(grep("Please run the model first, the trajectory is empty",
                      res[[1]]$message)) > 0)

res <- tools::assertError(trajectory(model, compartments = "I", as.is = TRUE))
stopifnot(length(grep("Please run the model first, the trajectory is empty",
                      res[[1]]$message)) > 0)

result <- run(model, threads = 1)
result

res <- tools::assertError(trajectory(result, I ~ S, as.is = TRUE))
stopifnot(res[[1]]$message == "Invalid formula specification of 'compartments'")

stopifnot(identical(length(trajectory(result, compartments = "S", as.is = TRUE)), 1001L))
stopifnot(identical(length(trajectory(result, compartments = "I", as.is = TRUE)), 1001L))

p <- prevalence(result, I~S+I, as.is = TRUE)
stopifnot(identical(length(p), 1001L))
stopifnot(is.null(dim(p)))

p <- prevalence(result, I~S+I, type = "wnp", as.is = TRUE)
stopifnot(identical(dim(p), c(1L, 1001L)))

if (SimInf:::have_openmp()) {
    result <- run(model, threads = 2)
    result

    stopifnot(identical(length(trajectory(result, compartments = "S", as.is = TRUE)), 1001L))
    stopifnot(identical(length(trajectory(result, compartments = "I", as.is = TRUE)), 1001L))

    p <- prevalence(result, I~S+I, as.is = TRUE)
    stopifnot(identical(length(p), 1001L))
    stopifnot(is.null(dim(p)))

    p <- prevalence(result, I~S+I, type = "wnp", as.is = TRUE)
    stopifnot(identical(dim(p), c(1L, 1001L)))
}

## Check measures for a SISe_sp model
model <- SISe_sp(u0       = data.frame(S = 99, I = 1),
                 tspan    = 0:1000,
                 events   = NULL,
                 phi      = 1,
                 upsilon  = 1,
                 gamma    = 0.1,
                 alpha    = 1,
                 beta_t1  = 1,
                 beta_t2  = 1,
                 beta_t3  = 1,
                 beta_t4  = 1,
                 end_t1   = 91,
                 end_t2   = 182,
                 end_t3   = 273,
                 end_t4   = 365,
                 coupling = 0,
                 distance = distance_matrix(1, 1, 1))

res <- tools::assertError(trajectory(model, compartments = "S", as.is = TRUE))
stopifnot(length(grep("Please run the model first, the trajectory is empty",
                      res[[1]]$message)) > 0)

res <- tools::assertError(trajectory(model, compartments = "I", as.is = TRUE))
stopifnot(length(grep("Please run the model first, the trajectory is empty",
                      res[[1]]$message)) > 0)

result <- run(model, threads = 1)
result

stopifnot(identical(length(trajectory(result, compartments = "S", as.is = TRUE)), 1001L))
stopifnot(identical(length(trajectory(result, compartments = "I", as.is = TRUE)), 1001L))

p <- prevalence(result, I~S+I, as.is = TRUE)
stopifnot(identical(length(p), 1001L))
stopifnot(is.null(dim(p)))

p <- prevalence(result, I~S+I, type = "wnp", as.is = TRUE)
stopifnot(identical(dim(p), c(1L, 1001L)))

if (SimInf:::have_openmp()) {
    result <- run(model, threads = 2)
    result

    stopifnot(identical(length(trajectory(result, compartments = "S", as.is = TRUE)), 1001L))
    stopifnot(identical(length(trajectory(result, compartments = "I", as.is = TRUE)), 1001L))

    p <- prevalence(result, I~S+I, as.is = TRUE)
    stopifnot(identical(length(p), 1001L))
    stopifnot(is.null(dim(p)))

    p <- prevalence(result, I~S+I, type = "wnp", as.is = TRUE)
    stopifnot(identical(dim(p), c(1L, 1001L)))
}

## Check 'susceptible' and 'infected' methods for a SISe3 model
u0 <- data.frame(S_1 = rep(10, 10), I_1 = rep( 0, 10),
                 S_2 = rep(20, 10), I_2 = rep( 0, 10),
                 S_3 = rep(70, 10), I_3 = rep( 0, 10))
model <- SISe3(u0        = u0,
               tspan     = seq_len(1000) - 1,
               events    = NULL,
               phi       = rep(1, 10),
               upsilon_1 = 0.0357,
               upsilon_2 = 0.0357,
               upsilon_3 = 0.00935,
               gamma_1   = 0.1,
               gamma_2   = 0.1,
               gamma_3   = 0.1,
               alpha     = 1.0,
               beta_t1   = 0.19,
               beta_t2   = 0.085,
               beta_t3   = 0.075,
               beta_t4   = 0.185,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 0.000011)

res <- tools::assertError(trajectory(model, compartments = "S_1", as.is = TRUE))
stopifnot(length(grep("Please run the model first, the trajectory is empty",
                      res[[1]]$message)) > 0)

res <- tools::assertError(trajectory(model, compartments = "I_1", as.is = TRUE))
stopifnot(length(grep("Please run the model first, the trajectory is empty",
                      res[[1]]$message)) > 0)

result <- run(model, threads = 1)
result

stopifnot(identical(length(trajectory(result, compartments = "S_1", as.is = TRUE)), 10000L))
stopifnot(identical(length(trajectory(result, compartments = "I_1", as.is = TRUE)), 10000L))

p <- prevalence(result, I_1+I_2+I_3~., as.is = TRUE)
stopifnot(identical(length(p), 1000L))
stopifnot(is.null(dim(p)))

p <- prevalence(result, I_1+I_2+I_3~., type = "wnp", as.is = TRUE)
stopifnot(identical(dim(p), c(10L, 1000L)))

if (SimInf:::have_openmp()) {
    result <- run(model, threads = 2)
    result

    stopifnot(identical(length(trajectory(result, compartments = "S_1", as.is = TRUE)), 10000L))
    stopifnot(identical(length(trajectory(result, compartments = "I_1", as.is = TRUE)), 10000L))

    p <- prevalence(result, I_1+I_2+I_3~., as.is = TRUE)
    stopifnot(identical(length(p), 1000L))
    stopifnot(is.null(dim(p)))

    p <- prevalence(result, I_1+I_2+I_3~., type = "wnp", as.is = TRUE)
    stopifnot(identical(dim(p), c(10L, 1000L)))
}

## Check measures with a SISe3_sp model
u0 <- data.frame(S_1 = rep(10, 10), I_1 = rep( 0, 10),
                 S_2 = rep(20, 10), I_2 = rep( 0, 10),
                 S_3 = rep(70, 10), I_3 = rep( 0, 10))

model <- SISe3_sp(u0        = u0,
                  tspan     = seq_len(1000) - 1,
                  events    = NULL,
                  phi       = rep(1, 10),
                  upsilon_1 = 0.0357,
                  upsilon_2 = 0.0357,
                  upsilon_3 = 0.00935,
                  gamma_1   = 0.1,
                  gamma_2   = 0.1,
                  gamma_3   = 0.1,
                  alpha     = 1.0,
                  beta_t1   = 0.19,
                  beta_t2   = 0.085,
                  beta_t3   = 0.075,
                  beta_t4   = 0.185,
                  end_t1    = 91,
                  end_t2    = 182,
                  end_t3    = 273,
                  end_t4    = 365,
                  coupling = 0,
                  distance = distance_matrix(1:10, 1:10, 1))

res <- tools::assertError(trajectory(model, compartments = "S_1", as.is = TRUE))
stopifnot(length(grep("Please run the model first, the trajectory is empty",
                      res[[1]]$message)) > 0)

res <- tools::assertError(trajectory(model, compartments = "I_1", as.is = TRUE))
stopifnot(length(grep("Please run the model first, the trajectory is empty",
                      res[[1]]$message)) > 0)

result <- run(model, threads = 1)
result

stopifnot(identical(length(trajectory(result, compartments = "S_1", as.is = TRUE)), 10000L))
stopifnot(identical(length(trajectory(result, compartments = "I_1", as.is = TRUE)), 10000L))

p <- prevalence(result, I_1+I_2+I_3~., as.is = TRUE)
stopifnot(identical(length(p), 1000L))
stopifnot(is.null(dim(p)))

p <- prevalence(result, I_1+I_2+I_3~., type = "wnp", as.is = TRUE)
stopifnot(identical(dim(p), c(10L, 1000L)))

if (SimInf:::have_openmp()) {
    result <- run(model, threads = 2)
    result

    stopifnot(identical(length(trajectory(result, compartments = "S_1", as.is = TRUE)), 10000L))
    stopifnot(identical(length(trajectory(result, compartments = "I_1", as.is = TRUE)), 10000L))

    p <- prevalence(result, I_1+I_2+I_3~., as.is = TRUE)
    stopifnot(identical(length(p), 1000L))
    stopifnot(is.null(dim(p)))

    p <- prevalence(result, I_1+I_2+I_3~., type = "wnp", as.is = TRUE)
    stopifnot(identical(dim(p), c(10L, 1000L)))
}
