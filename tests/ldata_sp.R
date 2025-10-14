## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2025 Stefan Widgren
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
library(Matrix)
library(tools)
source("util/check.R")

## For debugging
sessionInfo()

## Define a tolerance
tol <- 1e-8

## Local model parameters
l <- matrix(c(rep(91, 10), rep(182, 10), rep(273, 10), rep(365, 10)),
            nrow  = 4,
            byrow = TRUE)
storage.mode(l) <- "double"

## Distance matrix
d <- new("dgCMatrix",
         i = c(1L, 2L, 0L, 2L, 3L, 0L, 1L, 3L, 4L, 1L, 2L, 4L, 5L,
               2L, 3L, 5L, 6L, 3L, 4L, 6L, 7L, 4L, 5L, 7L, 8L, 5L,
               6L, 8L, 9L, 6L, 7L, 9L, 7L, 8L),
         p = c(0L, 2L, 5L, 9L, 13L, 17L, 21L, 25L, 29L, 32L, 34L),
         Dim = c(10L, 10L),
         Dimnames = list(NULL, NULL),
         x = c(1.4142135623731, 2.82842712474619, 1.4142135623731,
               1.4142135623731, 2.82842712474619, 2.82842712474619,
               1.4142135623731, 1.4142135623731, 2.82842712474619,
               2.82842712474619, 1.4142135623731, 1.4142135623731,
               2.82842712474619, 2.82842712474619, 1.4142135623731,
               1.4142135623731, 2.82842712474619, 2.82842712474619,
               1.4142135623731, 1.4142135623731, 2.82842712474619,
               2.82842712474619, 1.4142135623731, 1.4142135623731,
               2.82842712474619, 2.82842712474619, 1.4142135623731,
               1.4142135623731, 2.82842712474619, 2.82842712474619,
               1.4142135623731, 1.4142135623731, 2.82842712474619,
               1.4142135623731),
         factors = list())

## Check 'distance_matrix' method
d_obs <- distance_matrix(1:10, 1:10, 3)
stopifnot(is(d_obs, "dgCMatrix"))
stopifnot(identical(d_obs@i, d@i))
stopifnot(identical(d_obs@p, d@p))
stopifnot(all(abs(d_obs@x - d@x) < tol))

res <- assertError(distance_matrix(rep(1, 10), rep(1, 10), 3, "min_dist"))
check_error(res, "Invalid 'min_dist' argument. Please provide 'min_dist' > 0.")

res <- assertError(distance_matrix(rep(1, 10), rep(1, 10), 3, -1))
check_error(res, "Invalid 'min_dist' argument. Please provide 'min_dist' > 0.")

res <- assertError(distance_matrix(x = numeric(0), y = 1, cutoff = 2))
check_error(res, "'x' must be a numeric vector with length > 0.")

res <- assertError(distance_matrix(x = 1:3, y = 1:2, cutoff = 2))
check_error(res, "'y' must be a numeric vector with length 3.")

res <- assertError(distance_matrix(x = 1:3, y = 1:3, cutoff = -2))
check_error(res, "'cutoff' must be > 0.")

res <- assertError(distance_matrix(x = 1:3, y = 1:3, cutoff = Inf))
check_error(res, "'cutoff' must be > 0.")

res <- assertError(distance_matrix(x = 1:3, y = c(4, NA, 6), cutoff = 1))
check_error(res, "Invalid distance for i=0 and j=1.")

d_exp <- new("dgCMatrix",
             i = c(2L, 0L),
             p = c(0L, 1L, 1L, 2L),
             Dim = c(3L, 3L),
             Dimnames = list(NULL, NULL),
             x = c(2.828427125, 2.828427125),
             factors = list())
d_obs <- distance_matrix(x = 1:3, y = c(4, NA, 6), cutoff = 3, na_fail = FALSE)
stopifnot(is(d_obs, "dgCMatrix"))
stopifnot(identical(d_obs@i, d_exp@i))
stopifnot(identical(d_obs@p, d_exp@p))
stopifnot(all(abs(d_obs@x - d_exp@x) < tol))

res <- assertError(.Call(SimInf:::SimInf_distance_matrix,
                         x = c(1, 2, 3),
                         y = c(4, NA, 6),
                         cutoff = 3,
                         as.numeric(NULL),
                         na_fail = 1))
check_error(res, "'na_fail' must be TRUE or FALSE.")

## Check 'data' argument to C function 'SimInf_ldata_sp'
res <- assertError(.Call(SimInf:::SimInf_ldata_sp, NULL, d, 0L))
check_error(res, "Invalid 'data' argument.")

res <- assertError(.Call(SimInf:::SimInf_ldata_sp, d, d, 0L))
check_error(res, "Invalid 'data' argument.")

res <- assertError(.Call(SimInf:::SimInf_ldata_sp, 1:10, d, 0L))
check_error(res, "Invalid 'data' argument.")

## Check 'distance' argument to C function 'SimInf_ldata_sp'
res <- assertError(.Call(SimInf:::SimInf_ldata_sp, l, NULL, 0L))
check_error(res, "Invalid 'distance' argument.")

res <- assertError(.Call(SimInf:::SimInf_ldata_sp, l, l, 0L))
check_error(res, "Invalid 'distance' argument.")

res <- assertError(.Call(SimInf:::SimInf_ldata_sp, l, Diagonal(10), 0L))
check_error(res, "Invalid 'distance' argument.")

## Check 'metric' argument to C function 'SimInf_ldata_sp'
res <- assertError(.Call(SimInf:::SimInf_ldata_sp, l, d, NA_integer_))
check_error(res, "Invalid 'metric' argument.")

res <- assertError(.Call(SimInf:::SimInf_ldata_sp, l, d, NULL))
check_error(res, "Invalid 'metric' argument.")

res <- assertError(.Call(SimInf:::SimInf_ldata_sp, l, d, 0.0))
check_error(res, "Invalid 'metric' argument.")

res <- assertError(.Call(SimInf:::SimInf_ldata_sp, l, d, c(0L, 0L)))
check_error(res, "Invalid 'metric' argument.")

## Check non-equal number of nodes in 'distance' and 'data'
res <- assertError(.Call(SimInf:::SimInf_ldata_sp, l[, -1], d, 0L))
check_error(res, "The number of nodes in 'data' and 'distance' are not equal.")

## Check 'ldata' with metric equal to degree
ldata_exp <- structure(c(91, 182, 273, 365, 1, 3, 2, 4, -1, 0, 0, 0, 0, 0,
                         91, 182, 273, 365, 0, 2, 2, 4, 3, 4, -1, 0, 0, 0,
                         91, 182, 273, 365, 0, 2, 1, 3, 3, 4, 4, 4, -1, 0,
                         91, 182, 273, 365, 1, 3, 2, 4, 4, 4, 5, 4, -1, 0,
                         91, 182, 273, 365, 2, 4, 3, 4, 5, 4, 6, 4, -1, 0,
                         91, 182, 273, 365, 3, 4, 4, 4, 6, 4, 7, 4, -1, 0,
                         91, 182, 273, 365, 4, 4, 5, 4, 7, 4, 8, 3, -1, 0,
                         91, 182, 273, 365, 5, 4, 6, 4, 8, 3, 9, 2, -1, 0,
                         91, 182, 273, 365, 6, 4, 7, 4, 9, 2, -1, 0, 0, 0,
                         91, 182, 273, 365, 7, 4, 8, 3, -1, 0, 0, 0, 0, 0),
                       .Dim = c(14L, 10L))
ldata_obs <- .Call(SimInf:::SimInf_ldata_sp, l, d, 0L)
stopifnot(all(abs(ldata_obs - ldata_exp) < tol))

## Check 'ldata' with metric equal to distance
ldata_exp <- structure(c(91, 182, 273, 365, 1, 1.4142135623731,
                         2, 2.82842712474619, -1, 0, 0, 0, 0, 0,
                         91, 182, 273, 365, 0, 1.4142135623731, 2,
                         1.4142135623731, 3, 2.82842712474619, -1, 0, 0, 0,
                         91, 182, 273, 365, 0, 2.82842712474619, 1,
                         1.4142135623731, 3, 1.4142135623731,
                         4, 2.82842712474619, -1, 0,
                         91, 182, 273, 365, 1, 2.82842712474619,
                         2, 1.4142135623731, 4, 1.4142135623731,
                         5, 2.82842712474619, -1, 0,
                         91, 182, 273, 365, 2, 2.82842712474619,
                         3, 1.4142135623731, 5, 1.4142135623731,
                         6, 2.82842712474619, -1, 0,
                         91, 182, 273, 365, 3, 2.82842712474619,
                         4, 1.4142135623731, 6, 1.4142135623731,
                         7, 2.82842712474619, -1, 0,
                         91, 182, 273, 365, 4, 2.82842712474619,
                         5, 1.4142135623731, 7, 1.4142135623731,
                         8, 2.82842712474619, -1, 0,
                         91, 182, 273, 365, 5, 2.82842712474619,
                         6, 1.4142135623731, 8, 1.4142135623731,
                         9, 2.82842712474619, -1, 0,
                         91, 182, 273, 365, 6, 2.82842712474619,
                         7, 1.4142135623731, 9, 1.4142135623731,
                         -1, 0, 0, 0,
                         91, 182, 273, 365, 7, 2.82842712474619,
                         8, 1.4142135623731, -1, 0, 0, 0, 0, 0),
                       .Dim = c(14L, 10L))
ldata_obs <- add_spatial_coupling_to_ldata(x = 1:10, y = 1:10,
                                           cutoff = 3, ldata = l)
stopifnot(all(abs(ldata_obs - ldata_exp) < tol))

ldata_obs <- add_spatial_coupling_to_ldata(x = 1:10, y = 1:10, cutoff = 3)
stopifnot(all(abs(ldata_obs - ldata_exp[-(1:4), ]) < tol))

res <- assertError(add_spatial_coupling_to_ldata(x = 1:10,
                                                 y = 1:10,
                                                 cutoff = 3,
                                                 ldata = cbind(l, l)))
check_error(res, "Number of nodes in 'ldata' and coordinates must match.")

## Check 'ldata' with metric equal to 1 / distance^2
ldata_exp <- structure(c(91, 182, 273, 365, 1, 0.499999999999996, 2, 0.125,
                         -1, 0, 0, 0, 0, 0, 91, 182, 273, 365, 0,
                         0.499999999999996, 2, 0.499999999999996, 3, 0.125,
                         -1, 0, 0, 0, 91, 182, 273, 365, 0, 0.125, 1,
                         0.499999999999996, 3, 0.499999999999996, 4, 0.125,
                         -1, 0, 91, 182, 273, 365, 1, 0.125, 2,
                         0.499999999999996, 4, 0.499999999999996, 5, 0.125,
                         -1, 0, 91, 182, 273, 365, 2, 0.125, 3,
                         0.499999999999996, 5, 0.499999999999996, 6, 0.125,
                         -1, 0, 91, 182, 273, 365, 3, 0.125, 4,
                         0.499999999999996, 6, 0.499999999999996, 7, 0.125,
                         -1, 0, 91, 182, 273, 365, 4, 0.125, 5,
                         0.499999999999996, 7, 0.499999999999996, 8, 0.125,
                         -1, 0, 91, 182, 273, 365, 5, 0.125, 6,
                         0.499999999999996, 8, 0.499999999999996, 9, 0.125,
                         -1, 0, 91, 182, 273, 365, 6, 0.125, 7,
                         0.499999999999996, 9, 0.499999999999996, -1, 0, 0,
                         0, 91, 182, 273, 365, 7, 0.125, 8, 0.499999999999996,
                         -1, 0, 0, 0, 0, 0), .Dim = c(14L, 10L))

ldata_obs <- .Call(SimInf:::SimInf_ldata_sp, l, d, 2L)
stopifnot(all(abs(ldata_obs - ldata_exp) < tol))

## Check identical coordinates
res <- assertError(
    distance_matrix(x = c(1, 10, 1), y = c(1, 10, 1), cutoff = 20))
check_error(res, "Invalid 'min_dist' argument. Please provide 'min_dist' > 0.")

d_exp <- new("dgCMatrix",
             i = c(1L, 2L, 0L, 2L, 0L, 1L),
             p = c(0L, 2L, 4L, 6L),
             Dim = c(3L, 3L),
             Dimnames = list(NULL, NULL),
             x = c(12.7279220613579, 2, 12.7279220613579,
                   12.7279220613579, 2, 12.7279220613579),
         factors = list())
d_obs <- distance_matrix(x = c(1, 10, 1), y = c(1, 10, 1),
                         cutoff = 20, min_dist = 2)
stopifnot(is(d_obs, "dgCMatrix"))
stopifnot(identical(d_obs@i, d_exp@i))
stopifnot(identical(d_obs@p, d_exp@p))
stopifnot(all(abs(d_obs@x - d_exp@x) < tol))
