## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
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

## For debugging
sessionInfo()

## 1 dense discrete matrix.
## 0 dense continuous matrix.
## 1 replicate.
df_expected <- data.frame(
    node = rep(1:6, 10),
    time = rep(1:10, each = 6),
    S = seq(from = 1L, by = 3L, length.out = 60L),
    I = seq(from = 2L, by = 3L, length.out = 60L),
    R = seq(from = 3L, by = 3L, length.out = 60L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(1:180, nrow = 18, ncol = 10),
    1:3,
    c("S", "I", "R"),
    matrix(numeric(0), nrow = 0, ncol = 10),
    integer(0),
    NULL,
    as.numeric(1:10),
    6L,
    NULL,
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## 1 dense discrete matrix. Subset of columns.
## 0 dense continuous matrix.
## 1 replicate.
df_expected <- data.frame(
    node = rep(1:6, 10),
    time = rep(1:10, each = 6),
    S = seq(from = 1L, by = 3L, length.out = 60L),
    R = seq(from = 3L, by = 3L, length.out = 60L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(1:180, nrow = 18, ncol = 10),
    c(1L, 3L),
    c("S", "I", "R"),
    matrix(numeric(0), nrow = 0, ncol = 10),
    integer(0),
    NULL,
    as.numeric(1:10),
    6L,
    NULL,
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## 1 dense discrete matrix. Subset of columns. Subset of identifiers.
## 0 dense continuous matrix.
## 1 replicate.
df_expected <- data.frame(
    node = rep(c(1L, 3L, 6L), 10),
    time = rep(1:10, each = 3),
    S = c(1L, 7L, 16L, 19L, 25L, 34L, 37L, 43L, 52L, 55L, 61L,
          70L, 73L, 79L, 88L, 91L, 97L, 106L, 109L, 115L, 124L,
          127L, 133L, 142L, 145L, 151L, 160L, 163L, 169L, 178L),
    R = c(3L, 9L, 18L, 21L, 27L, 36L, 39L, 45L, 54L, 57L, 63L,
          72L, 75L, 81L, 90L, 93L, 99L, 108L, 111L, 117L, 126L,
          129L, 135L, 144L, 147L, 153L, 162L, 165L, 171L, 180L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(1:180, nrow = 18, ncol = 10),
    c(1L, 3L),
    c("S", "I", "R"),
    matrix(numeric(0), nrow = 0, ncol = 10),
    integer(0),
    NULL,
    as.numeric(1:10),
    6L,
    c(1L, 3L, 6L),
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## 1 dense discrete matrix.
## 1 dense continuous matrix.
## 1 replicate.
df_expected <- data.frame(
    node = rep(1:6, 10),
    time = rep(1:10, each = 6),
    S = seq(from = 1L, by = 3L, length.out = 60L),
    I = seq(from = 2L, by = 3L, length.out = 60L),
    R = seq(from = 3L, by = 3L, length.out = 60L),
    phi = seq(from = 1000, by = 1, length.out = 60L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(1:180, nrow = 18, ncol = 10),
    1:3,
    c("S", "I", "R"),
    matrix(as.numeric(1000:1059), nrow = 6, ncol = 10),
    1L,
    "phi",
    as.numeric(1:10),
    6L,
    NULL,
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## 0 dense discrete matrix.
## 1 dense continuous matrix.
## 1 replicate.
df_expected <- data.frame(
    node = rep(1:6, 10),
    time = rep(1:10, each = 6),
    phi = seq(from = 1000, by = 1, length.out = 60L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(integer(0), nrow = 0, ncol = 10),
    integer(0),
    NULL,
    matrix(as.numeric(1000:1059), nrow = 6, ncol = 10),
    1L,
    "phi",
    as.numeric(1:10),
    6L,
    NULL,
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## 1 dense discrete matrix.
## 0 dense continuous matrix.
## 2 replicates.
df_expected <- data.frame(
    node = rep(rep(1:6, 10), 2),
    time = rep(rep(1:10, each = 6), 2),
    replicate = rep(1:2, each = 60),
    S = seq(from = 1L, by = 3L, length.out = 120L),
    I = seq(from = 2L, by = 3L, length.out = 120L),
    R = seq(from = 3L, by = 3L, length.out = 120L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(1:360, nrow = 18, ncol = 20),
    1:3,
    c("S", "I", "R"),
    matrix(numeric(0), nrow = 0, ncol = 20),
    integer(0),
    NULL,
    as.numeric(1:10),
    6L,
    NULL,
    "node",
    2L)
stopifnot(identical(df_observed, df_expected))
