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
library(Matrix)

## For debugging
sessionInfo()

## 1 dense discrete matrix.
## No continuous matrix.
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

## 1 dense discrete matrix.
## No continuous matrix.
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

## 1 dense discrete matrix.
## No continuous matrix.
## 1 replicate.
df_expected <- data.frame(
    node = rep(1:6, 10),
    time = rep(c("2025-01-01", "2025-01-02", "2025-01-03", "2025-01-04",
                 "2025-01-05", "2025-01-06", "2025-01-07", "2025-01-08",
                 "2025-01-09", "2025-01-10"),
               each = 6),
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
    c("2025-01-01" = 1, "2025-01-02" = 2, "2025-01-03" = 3,
      "2025-01-04" = 4, "2025-01-05" = 5, "2025-01-06" = 6,
      "2025-01-07" = 7, "2025-01-08" = 8, "2025-01-09" = 9,
      "2025-01-10" = 10),
    6L,
    NULL,
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## 1 dense discrete matrix. Subset of columns.
## No continuous matrix.
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

## 1 dense discrete matrix. Subset of columns.
## No continuous matrix.
## 2 replicates.
df_expected <- data.frame(
    node = rep(rep(1:6, 10), 2),
    time = rep(rep(1:10, each = 6), 2),
    replicate = rep(1:2, each = 60),
    S = seq(from = 1L, by = 3L, length.out = 120L),
    R = seq(from = 3L, by = 3L, length.out = 120L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(1:360, nrow = 18, ncol = 20),
    c(1L, 3L),
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

## 1 dense discrete matrix. Subset of columns. Subset of identifiers.
## No continuous matrix.
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

## 1 dense discrete matrix. Subset of columns. Subset of identifiers.
## No continuous matrix.
## 2 replicates.
df_expected <- data.frame(
    node = rep(rep(c(1L, 3L, 6L), 10), 2),
    time = rep(rep(1:10, each = 3), 2),
    replicate = rep(1:2, each = 30),
    S = c(1L, 7L, 16L, 19L, 25L, 34L, 37L, 43L, 52L, 55L, 61L,
          70L, 73L, 79L, 88L, 91L, 97L, 106L, 109L, 115L, 124L,
          127L, 133L, 142L, 145L, 151L, 160L, 163L, 169L, 178L,
          181L, 187L, 196L, 199L, 205L, 214L, 217L, 223L, 232L,
          235L, 241L, 250L, 253L, 259L, 268L, 271L, 277L, 286L,
          289L, 295L, 304L, 307L, 313L, 322L, 325L, 331L, 340L,
          343L, 349L, 358L),
    R = c(3L, 9L, 18L, 21L, 27L, 36L, 39L, 45L, 54L, 57L, 63L,
          72L, 75L, 81L, 90L, 93L, 99L, 108L, 111L, 117L, 126L,
          129L, 135L, 144L, 147L, 153L, 162L, 165L, 171L, 180L,
          183L, 189L, 198L, 201L, 207L, 216L, 219L, 225L, 234L,
          237L, 243L, 252L, 255L, 261L, 270L, 273L, 279L, 288L,
          291L, 297L, 306L, 309L, 315L, 324L, 327L, 333L, 342L,
          345L, 351L, 360L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(1:360, nrow = 18, ncol = 20),
    c(1L, 3L),
    c("S", "I", "R"),
    matrix(numeric(0), nrow = 0, ncol = 20),
    integer(0),
    NULL,
    as.numeric(1:10),
    6L,
    c(1L, 3L, 6L),
    "node",
    2L)
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

## 1 dense discrete matrix.
## 1 dense continuous matrix.
## 2 replicates.
df_expected <- data.frame(
    node = rep(rep(1:6, 10), 2),
    time = rep(rep(1:10, each = 6), 2),
    replicate = rep(1:2, each = 60),
    S = seq(from = 1L, by = 3L, length.out = 120L),
    I = seq(from = 2L, by = 3L, length.out = 120L),
    R = seq(from = 3L, by = 3L, length.out = 120L),
    phi = seq(from = 1000, by = 1, length.out = 120L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(1:360, nrow = 18, ncol = 20),
    1:3,
    c("S", "I", "R"),
    matrix(as.numeric(1000:1119), nrow = 6, ncol = 20),
    1L,
    "phi",
    as.numeric(1:10),
    6L,
    NULL,
    "node",
    2L)
stopifnot(identical(df_observed, df_expected))

## No discrete matrix.
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

## No discrete matrix.
## 1 dense continuous matrix.
## 2 replicates.
df_expected <- data.frame(
    node = rep(rep(1:6, 10), 2),
    time = rep(rep(1:10, each = 6), 2),
    replicate = rep(1:2, each = 60),
    phi = seq(from = 1000, by = 1, length.out = 120L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(integer(0), nrow = 0, ncol = 20),
    integer(0),
    NULL,
    matrix(as.numeric(1000:1119), nrow = 6, ncol = 20),
    1L,
    "phi",
    as.numeric(1:10),
    6L,
    NULL,
    "node",
    2L)
stopifnot(identical(df_observed, df_expected))

## 1 sparse discrete matrix.
## No continuous matrix.
## 1 replicate.
dm <- matrix(1:180, nrow = 18, ncol = 10)
dm[c(4, 5, 6, 13, 14, 15), ] <- 0L
dm[1, 1] <- 0L
dm[2, 2] <- 0L
dm[3, 3] <- 0L
dm[7, 4] <- 0L
dm[8, 5] <- 0L
dm[9, 6] <- 0L
dm[10, 7] <- 0L
dm[11, 8] <- 0L
dm[12, 9] <- 0L
dm[16, 10] <- 0L
dm[17, 1] <- 0L
dm[18, 1] <- 0L
dm <- as(as(as(Matrix(dm), "dMatrix"), "generalMatrix"), "CsparseMatrix")

df_expected <- data.frame(
    node = c(1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L,
             1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L,
             1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L),
    time = c("2025-01-01", "2025-01-01", "2025-01-01", "2025-01-01",
             "2025-01-02", "2025-01-02", "2025-01-02", "2025-01-02",
             "2025-01-03", "2025-01-03", "2025-01-03", "2025-01-03",
             "2025-01-04", "2025-01-04", "2025-01-04", "2025-01-04",
             "2025-01-05", "2025-01-05", "2025-01-05", "2025-01-05",
             "2025-01-06", "2025-01-06", "2025-01-06", "2025-01-06",
             "2025-01-07", "2025-01-07", "2025-01-07", "2025-01-07",
             "2025-01-08", "2025-01-08", "2025-01-08", "2025-01-08",
             "2025-01-09", "2025-01-09", "2025-01-09", "2025-01-09",
             "2025-01-10", "2025-01-10", "2025-01-10", "2025-01-10"),
    S = c(NA, 7L, 10L, 16L, 19L, 25L, 28L, 34L, 37L, 43L, 46L, 52L, 55L, NA,
          64L, 70L, 73L, 79L, 82L, 88L, 91L, 97L, 100L, 106L, 109L, 115L,
          NA, 124L, 127L, 133L, 136L, 142L, 145L, 151L, 154L, 160L, 163L,
          169L, 172L, NA),
    I = c(2L, 8L, 11L, NA, NA, 26L, 29L, 35L, 38L, 44L, 47L, 53L, 56L, 62L,
          65L, 71L, 74L, NA, 83L, 89L, 92L, 98L, 101L, 107L, 110L, 116L,
          119L, 125L, 128L, 134L, NA, 143L, 146L, 152L, 155L, 161L, 164L,
          170L, 173L, 179L),
    R = c(3L, 9L, 12L, NA, 21L, 27L, 30L, 36L, NA, 45L, 48L, 54L, 57L, 63L,
          66L, 72L, 75L, 81L, 84L, 90L, 93L, NA, 102L, 108L, 111L, 117L,
          120L, 126L, 129L, 135L, 138L, 144L, 147L, 153L, NA, 162L, 165L,
          171L, 174L, 180L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    dm,
    1:3,
    c("S", "I", "R"),
    matrix(numeric(0), nrow = 0, ncol = 10),
    integer(0),
    NULL,
    c("2025-01-01" = 1, "2025-01-02" = 2, "2025-01-03" = 3,
      "2025-01-04" = 4, "2025-01-05" = 5, "2025-01-06" = 6,
      "2025-01-07" = 7, "2025-01-08" = 8, "2025-01-09" = 9,
      "2025-01-10" = 10),
    6L,
    NULL,
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## 1 sparse discrete matrix. Subset of columns.
## No continuous matrix.
## 1 replicate.

df_expected <- data.frame(
    node = c(1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L, 1L, 3L,
             4L, 6L, 1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L,
             1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L, 1L, 3L, 4L, 6L),
    time = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 4L, 4L,
             4L, 4L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 7L, 7L, 7L, 7L,
             8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 10L),
    S = c(NA, 7L, 10L, 16L, 19L, 25L, 28L, 34L, 37L, 43L, 46L, 52L,
          55L, NA, 64L, 70L, 73L, 79L, 82L, 88L, 91L, 97L, 100L,
          106L, 109L, 115L, NA, 124L, 127L, 133L, 136L, 142L, 145L,
          151L, 154L, 160L, 163L, 169L, 172L, NA),
    R = c(3L, 9L, 12L, NA, 21L, 27L, 30L, 36L, NA, 45L, 48L, 54L, 57L,
          63L, 66L, 72L, 75L, 81L, 84L, 90L, 93L, NA, 102L, 108L, 111L,
          117L, 120L, 126L, 129L, 135L, 138L, 144L, 147L, 153L, NA,
          162L, 165L, 171L, 174L, 180L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    dm,
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

## 1 sparse discrete matrix. Subset of identifiers.
## No continuous matrix.
## 1 replicate.

df_expected <- data.frame(
    node = c(1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L,
             6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L,
             3L, 6L),
    time = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L, 5L, 5L,
             5L, 6L, 6L, 6L, 7L, 7L, 7L, 8L, 8L, 8L, 9L, 9L, 9L, 10L,
             10L, 10L),
    S = c(NA, 7L, 16L, 19L, 25L, 34L, 37L, 43L, 52L, 55L, NA, 70L,
          73L, 79L, 88L, 91L, 97L, 106L, 109L, 115L, 124L, 127L,
          133L, 142L, 145L, 151L, 160L, 163L, 169L, NA),
    I = c(2L, 8L, NA, NA, 26L, 35L, 38L, 44L, 53L, 56L, 62L, 71L,
          74L, NA, 89L, 92L, 98L, 107L, 110L, 116L, 125L, 128L, 134L,
          143L, 146L, 152L, 161L, 164L, 170L, 179L),
    R = c(3L, 9L, NA, 21L, 27L, 36L, NA, 45L, 54L, 57L, 63L, 72L,
          75L, 81L, 90L, 93L, NA, 108L, 111L, 117L, 126L, 129L, 135L,
          144L, 147L, 153L, 162L, 165L, 171L, 180L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    dm,
    1:3,
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

## 1 sparse discrete matrix. Subset of columns. Subset of identifiers.
## No continuous matrix.
## 1 replicate.

df_expected <- data.frame(
    node = c(1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L,
             6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L,
             3L, 6L),
    time = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L, 5L, 5L,
             5L, 6L, 6L, 6L, 7L, 7L, 7L, 8L, 8L, 8L, 9L, 9L, 9L, 10L,
             10L, 10L),
    S = c(NA, 7L, 16L, 19L, 25L, 34L, 37L, 43L, 52L, 55L, NA, 70L,
          73L, 79L, 88L, 91L, 97L, 106L, 109L, 115L, 124L, 127L, 133L,
          142L, 145L, 151L, 160L, 163L, 169L, NA),
    R = c(3L, 9L, NA, 21L, 27L, 36L, NA, 45L, 54L, 57L, 63L, 72L, 75L,
          81L, 90L, 93L, NA, 108L, 111L, 117L, 126L, 129L, 135L, 144L,
          147L, 153L, 162L, 165L, 171L, 180L))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    dm,
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

## No discrete matrix.
## 1 sparse continuous matrix.
## 1 replicate.
cm <- matrix(as.numeric(1001:1180), nrow = 18, ncol = 10)
cm[c(13, 14, 15), ] <- 0L
cm[1, 1] <- 0L
cm[2, 2] <- 0L
cm[3, 3] <- 0L
cm[7, 4] <- 0L
cm[8, 5] <- 0L
cm[9, 6] <- 0L
cm[10, 7] <- 0L
cm[11, 8] <- 0L
cm[12, 9] <- 0L
cm[16, 10] <- 0L
cm[17, 1] <- 0L
cm[18, 1] <- 0L
cm <- as(as(as(Matrix(cm), "dMatrix"), "generalMatrix"), "CsparseMatrix")

df_expected <- data.frame(
    node = c(1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L,
             1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L,
             1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L,
             1L, 2L, 3L, 4L, 6L),
    time = c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L,
             4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L,
             7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L,
             10L, 10L, 10L, 10L, 10L),
    alpha = c(NA, 1004, 1007, 1010, 1016, 1019, 1022, 1025, 1028, 1034,
              1037, 1040, 1043, 1046, 1052, 1055, 1058, NA, 1064, 1070,
              1073, 1076, 1079, 1082, 1088, 1091, 1094, 1097, 1100, 1106,
              1109, 1112, 1115, NA, 1124, 1127, 1130, 1133, 1136, 1142,
              1145, 1148, 1151, 1154, 1160, 1163, 1166, 1169, 1172, NA),
    beta = c(1002, 1005, 1008, 1011, NA, NA, 1023, 1026, 1029, 1035,
             1038, 1041, 1044, 1047, 1053, 1056, 1059, 1062, 1065, 1071,
             1074, 1077, NA, 1083, 1089, 1092, 1095, 1098, 1101, 1107,
             1110, 1113, 1116, 1119, 1125, 1128, 1131, 1134, NA, 1143,
             1146, 1149, 1152, 1155, 1161, 1164, 1167, 1170, 1173, 1179),
    gamma = c(1003, 1006, 1009, 1012, NA, 1021, 1024, 1027, 1030, 1036,
              NA, 1042, 1045, 1048, 1054, 1057, 1060, 1063, 1066, 1072,
              1075, 1078, 1081, 1084, 1090, 1093, 1096, NA, 1102, 1108,
              1111, 1114, 1117, 1120, 1126, 1129, 1132, 1135, 1138, 1144,
              1147, 1150, 1153, NA, 1162, 1165, 1168, 1171, 1174, 1180))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(integer(0), nrow = 0, ncol = 10),
    integer(0),
    NULL,
    cm,
    1:3,
    c("alpha", "beta", "gamma"),
    as.numeric(1:10),
    6L,
    NULL,
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## No discrete matrix.
## 1 sparse continuous matrix. Subset of columns.
## 1 replicate.

df_expected <- data.frame(
    node = c(1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L,
             1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L,
             1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L,
             1L, 2L, 3L, 4L, 6L),
    time = c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L,
             4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L,
             7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L,
             10L, 10L, 10L, 10L, 10L),
    alpha = c(NA, 1004, 1007, 1010, 1016, 1019, 1022, 1025, 1028, 1034,
              1037, 1040, 1043, 1046, 1052, 1055, 1058, NA, 1064, 1070,
              1073, 1076, 1079, 1082, 1088, 1091, 1094, 1097, 1100, 1106,
              1109, 1112, 1115, NA, 1124, 1127, 1130, 1133, 1136, 1142,
              1145, 1148, 1151, 1154, 1160, 1163, 1166, 1169, 1172, NA),
    gamma = c(1003, 1006, 1009, 1012, NA, 1021, 1024, 1027, 1030, 1036,
              NA, 1042, 1045, 1048, 1054, 1057, 1060, 1063, 1066, 1072,
              1075, 1078, 1081, 1084, 1090, 1093, 1096, NA, 1102, 1108,
              1111, 1114, 1117, 1120, 1126, 1129, 1132, 1135, 1138, 1144,
              1147, 1150, 1153, NA, 1162, 1165, 1168, 1171, 1174, 1180))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(integer(0), nrow = 0, ncol = 10),
    integer(0),
    NULL,
    cm,
    c(1L, 3L),
    c("alpha", "beta", "gamma"),
    as.numeric(1:10),
    6L,
    NULL,
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## No discrete matrix.
## 1 sparse continuous matrix. Subset of identifiers.
## 1 replicate.
df_expected <- data.frame(
    node = c(1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L,
             1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L,
             1L, 3L, 6L, 1L, 3L, 6L),
    time = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L,
             5L, 5L, 5L, 6L, 6L, 6L, 7L, 7L, 7L, 8L, 8L, 8L,
             9L, 9L, 9L, 10L, 10L, 10L),
    alpha = c(NA, 1007, 1016, 1019, 1025,  1034, 1037, 1043,
              1052, 1055, NA, 1070, 1073, 1079, 1088, 1091,
              1097, 1106, 1109, 1115, 1124, 1127, 1133, 1142,
              1145, 1151, 1160, 1163, 1169, NA),
    beta = c(1002, 1008, NA, NA, 1026, 1035, 1038, 1044, 1053,
             1056, 1062, 1071, 1074, NA, 1089, 1092, 1098, 1107,
             1110, 1116, 1125, 1128, 1134, 1143, 1146, 1152,
             1161, 1164, 1170, 1179),
    gamma = c(1003, 1009, NA, 1021, 1027, 1036, NA, 1045, 1054,
              1057, 1063, 1072, 1075, 1081, 1090, 1093, NA,
              1108, 1111, 1117, 1126, 1129, 1135, 1144, 1147,
              1153, 1162, 1165, 1171, 1180))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(integer(0), nrow = 0, ncol = 10),
    integer(0),
    NULL,
    cm,
    1:3,
    c("alpha", "beta", "gamma"),
    as.numeric(1:10),
    6L,
    c(1L, 3L, 6L),
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## No discrete matrix.
## 1 sparse continuous matrix. Subset of columns and identifiers.
## 1 replicate.
df_expected <- data.frame(
    node = c(1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L,
             1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L, 1L, 3L, 6L,
             1L, 3L, 6L, 1L, 3L, 6L),
    time = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L,  4L, 4L, 4L,
             5L, 5L, 5L, 6L, 6L, 6L, 7L, 7L, 7L, 8L, 8L, 8L,
             9L, 9L, 9L, 10L, 10L, 10L),
    alpha = c(NA, 1007, 1016, 1019, 1025, 1034, 1037, 1043,
              1052, 1055, NA, 1070, 1073, 1079, 1088, 1091,
              1097, 1106, 1109, 1115, 1124, 1127, 1133, 1142,
              1145, 1151, 1160, 1163, 1169, NA),
    gamma = c(1003, 1009, NA, 1021, 1027, 1036, NA, 1045, 1054,
              1057, 1063, 1072, 1075, 1081, 1090, 1093, NA, 1108,
              1111, 1117, 1126, 1129, 1135, 1144, 1147, 1153,
              1162, 1165, 1171, 1180))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    matrix(integer(0), nrow = 0, ncol = 10),
    integer(0),
    NULL,
    cm,
    c(1L, 3L),
    c("alpha", "beta", "gamma"),
    as.numeric(1:10),
    6L,
    c(1L, 3L, 6L),
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))

## 1 sparse discrete matrix.
## 1 sparse continuous matrix.
## 1 replicate.

df_expected <- data.frame(
    node = c(1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L,
             1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L,
             1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L,
             1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L,
             1L, 2L, 3L, 4L, 6L, 1L, 2L, 3L, 4L, 6L),
    time = c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L,
             3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L,
             5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L,
             7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L,
             9L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 10L, 10L),
    S = c(NA, NA, 7L, 10L, 16L, 19L, NA, 25L, 28L, 34L,
          37L, NA, 43L, 46L, 52L, 55L, NA, NA, 64L, 70L,
          73L, NA, 79L, 82L, 88L, 91L, NA, 97L, 100L,
          106L, 109L, NA, 115L, NA, 124L, 127L, NA, 133L,
          136L, 142L, 145L, NA, 151L, 154L, 160L, 163L,
          NA, 169L, 172L, NA),
    I = c(2L, NA, 8L, 11L, NA, NA, NA, 26L, 29L, 35L, 38L,
          NA, 44L, 47L, 53L, 56L, NA, 62L, 65L, 71L, 74L,
          NA, NA, 83L, 89L, 92L, NA, 98L, 101L, 107L, 110L,
          NA, 116L, 119L, 125L, 128L, NA, 134L, NA, 143L,
          146L, NA, 152L, 155L, 161L, 164L, NA, 170L, 173L, 179L),
    R = c(3L, NA, 9L, 12L, NA, 21L, NA, 27L, 30L, 36L, NA,
          NA, 45L, 48L, 54L, 57L, NA, 63L, 66L, 72L, 75L, NA,
          81L, 84L, 90L, 93L, NA, NA, 102L, 108L, 111L,
          NA, 117L, 120L, 126L, 129L, NA, 135L, 138L, 144L,
          147L, NA, 153L, NA, 162L, 165L, NA, 171L, 174L, 180L),
    alpha = c(NA, 1004, 1007, 1010, 1016, 1019, 1022, 1025, 1028,
              1034, 1037, 1040, 1043, 1046, 1052, 1055, 1058,
              NA, 1064, 1070, 1073, 1076, 1079, 1082, 1088, 1091,
              1094, 1097, 1100, 1106, 1109, 1112, 1115, NA, 1124,
              1127, 1130, 1133, 1136, 1142, 1145, 1148, 1151,
              1154, 1160, 1163, 1166, 1169, 1172, NA),
    beta = c(1002, 1005, 1008, 1011, NA, NA, 1023, 1026, 1029,
             1035, 1038, 1041, 1044, 1047, 1053, 1056, 1059,
             1062, 1065, 1071, 1074, 1077, NA, 1083, 1089, 1092,
             1095, 1098, 1101, 1107, 1110, 1113, 1116, 1119, 1125,
             1128, 1131, 1134, NA, 1143, 1146, 1149, 1152, 1155,
             1161, 1164, 1167, 1170, 1173, 1179),
    gamma = c(1003, 1006, 1009, 1012, NA, 1021, 1024, 1027, 1030,
              1036, NA, 1042, 1045, 1048, 1054, 1057, 1060, 1063,
              1066, 1072, 1075, 1078, 1081, 1084, 1090, 1093,
              1096, NA, 1102, 1108, 1111, 1114, 1117, 1120, 1126,
              1129, 1132, 1135, 1138, 1144, 1147, 1150, 1153, NA,
              1162, 1165, 1168, 1171, 1174, 1180))
df_observed <- .Call(
    SimInf:::SimInf_trajectory,
    dm,
    1:3,
    c("S", "I", "R"),
    cm,
    1:3,
    c("alpha", "beta", "gamma"),
    as.numeric(1:10),
    6L,
    NULL,
    "node",
    1L)
stopifnot(identical(df_observed, df_expected))
