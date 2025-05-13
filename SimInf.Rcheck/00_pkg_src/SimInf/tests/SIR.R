## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2021 Stefan Widgren
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
library(tools)
source("util/check.R")

## Specify the number of threads to use.
set_num_threads(1)

## For debugging
sessionInfo()

## Check invalid u0
res <- assertError(SIR(u0 = "u0"))
check_error(res, "Missing columns in u0.")

u0 <- data.frame(S  = c(0, 1, 2, 3, 4, 5),
                 I  = c(0, 0, 0, 0, 0, 0),
                 R  = c(0, 0, 0, 0, 0, 0))

## Check missing columns in u0
res <- assertError(SIR(u0 = u0[, c("I", "R"), drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- assertError(SIR(u0 = u0[, c("S", "R"), drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- assertError(SIR(u0 = u0[, c("S", "I"), drop = FALSE]))
check_error(res, "Missing columns in u0.")

## Check missing beta
res <- assertError(SIR(u0     = u0,
                       tspan  = seq_len(10) - 1,
                       events = NULL,
                       gamma  = 0.5))
check_error(res, "'beta' must be numeric of length 1 or 'nrow(u0)'.")

## Check missing gamma
res <- assertError(SIR(u0     = u0,
                       tspan  = seq_len(10) - 1,
                       events = NULL,
                       beta   = 0.5))
check_error(res, "'gamma' must be numeric of length 1 or 'nrow(u0)'.")

## Check non-numeric beta
res <- assertError(SIR(u0      = u0,
                       tspan   = seq_len(10) - 1,
                       events  = NULL,
                       beta    = "0.5",
                       gamma   = 0.1))
check_error(res, "'beta' must be numeric of length 1 or 'nrow(u0)'.")

## Check non-numeric gamma
res <- assertError(SIR(u0      = u0,
                       tspan   = seq_len(10) - 1,
                       events  = NULL,
                       beta    = 0.5,
                       gamma   = "0.1"))
check_error(res, "'gamma' must be numeric of length 1 or 'nrow(u0)'.")

## Check that length of beta equals 1 or nrow(u0)
res <- assertError(SIR(u0      = u0,
                       tspan   = seq_len(10) - 1,
                       events  = NULL,
                       beta    = c(0.5, 0.5),
                       gamma   = 0.1))
check_error(res, "'beta' must be numeric of length 1 or 'nrow(u0)'.")

## Check that length of gamma equals 1 or nrow(u0)
res <- assertError(SIR(u0      = u0,
                       tspan   = seq_len(10) - 1,
                       events  = NULL,
                       beta    = 0.5,
                       gamma   = c(0.1, 0.1)))
check_error(res, "'gamma' must be numeric of length 1 or 'nrow(u0)'.")

## Extract data from the 'suscpetible', 'infected' and 'recovered'
## compartments
model <- SIR(u0     = u0,
             tspan  = seq_len(10) - 1,
             events = NULL,
             beta   = 0,
             gamma  = 0)

result <- run(model)

S_expected <- structure(c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L,
                          1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
                          2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L,
                          3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
                          4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                        .Dim = c(6L, 10L))

S_observed <- trajectory(result, compartments = "S", format = "matrix")
stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L))

I_observed <- trajectory(result, compartments = "I", format = "matrix")
stopifnot(identical(I_observed, I_expected))

R_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L))

R_observed <- trajectory(result, compartments = "R", format = "matrix")

stopifnot(identical(R_observed, R_expected))

## Check data
stopifnot(identical(nrow(events_SIR()), 466692L))
stopifnot(identical(nrow(u0_SIR()), 1600L))

## Check SIR plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SIR prevalence plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result, I ~ S + I + R)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SIR events plot with no events
model <- SIR(u0     = u0,
             tspan  = seq_len(10) - 1,
             events = NULL,
             beta   = 0,
             gamma  = 0)
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(model@events)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SIR plot with tspan as Date vector
model <- SIR(u0     = u0,
             tspan  = as.Date(seq_len(10), origin = "2016-12-31"),
             events = NULL,
             beta   = 0,
             gamma  = 0)
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(run(model))
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SIR events plot method
model <- SIR(u0     = u0_SIR(),
             tspan  = seq_len(365 * 4),
             events = events_SIR(),
             beta   = 0,
             gamma  = 0)
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(model@events)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SIR run function fails for misspecified SIR model
res <- assertError(.Call(SimInf:::SIR_run, NULL, NULL))
check_error(res, "Invalid model.")

res <- assertError(.Call(SimInf:::SIR_run, "SIR", NULL))
check_error(res, "Invalid model.")

model <- SIR(u0     = u0_SIR(),
             tspan  = seq_len(365 * 4),
             events = events_SIR(),
             beta   = 0,
             gamma  = 0)
stopifnot(is(events(model), "SimInf_events"))

id <- c(64L, 63L, 62L, 61L, 64L, 57L, 61L, 65L, 70L, 71L, 52L, 56L,
        73L, 68L, 74L, 78L, 67L, 64L, 50L, 67L, 62L, 69L, 66L, 59L,
        60L, 55L, 55L, 61L, 59L, 56L, 59L, 56L, 77L, 71L, 54L, 63L,
        61L, 78L, 77L, 81L, 47L, 64L, 72L, 57L, 63L, 57L, 55L, 63L,
        48L, 82L, 62L, 70L, 67L, 62L, 68L, 54L, 59L, 68L, 60L, 57L,
        57L, 48L, 68L, 58L, 75L, 57L, 73L, 68L, 66L, 66L, 68L, 68L,
        63L, 63L, 60L, 66L, 64L, 78L, 70L, 57L, 68L, 55L, 63L, 66L,
        59L, 58L, 64L, 57L, 57L, 52L, 67L, 56L, 72L, 65L, 57L, 55L,
        67L, 68L, 72L, 65L, 57L, 62L, 45L, 74L, 65L, 60L, 66L, 72L,
        59L, 65L, 57L, 55L, 58L, 67L, 55L, 60L, 65L, 56L, 71L, 60L,
        55L, 64L, 48L, 65L, 70L, 67L, 61L, 58L, 66L, 69L, 72L, 61L,
        60L, 69L, 57L, 58L, 59L, 59L, 66L, 72L, 70L, 54L, 73L, 55L,
        62L, 61L, 68L, 64L, 72L, 54L, 47L, 65L, 58L, 52L, 60L, 68L,
        56L, 71L, 68L, 66L, 72L, 59L, 62L, 71L, 52L, 65L, 80L, 71L,
        56L, 52L, 60L, 70L, 75L, 55L, 70L, 73L, 49L, 66L, 59L, 60L,
        61L, 60L, 74L, 62L, 70L, 51L, 63L, 48L, 63L, 64L, 55L, 67L,
        58L, 87L, 75L, 65L, 57L, 58L, 66L, 61L, 58L, 68L, 60L, 67L,
        66L, 61L, 58L, 59L, 74L, 44L, 65L, 75L, 63L, 60L, 59L, 65L,
        63L, 52L, 59L, 74L, 67L, 63L, 61L, 55L, 64L, 66L, 65L, 63L,
        61L, 62L, 53L, 60L, 79L, 60L, 58L, 64L, 70L, 78L, 64L, 70L,
        57L, 60L, 62L, 71L, 68L, 56L, 63L, 55L, 54L, 71L, 61L, 78L,
        64L, 53L, 49L, 43L, 60L, 45L, 65L, 59L, 71L, 59L, 63L, 60L,
        76L, 42L, 71L, 48L, 72L, 48L, 54L, 60L, 64L, 57L, 61L, 76L,
        75L, 68L, 63L, 56L, 66L, 63L, 56L, 58L, 71L, 70L, 57L, 49L,
        62L, 63L, 58L, 61L, 51L, 61L, 69L, 72L, 61L, 53L, 65L, 61L,
        68L, 61L, 61L, 71L, 61L, 58L, 64L, 73L, 44L, 72L, 51L, 48L,
        67L, 86L, 54L, 68L, 69L, 57L, 60L, 54L, 74L, 59L, 50L, 72L,
        53L, 54L, 47L, 57L, 73L, 58L, 55L, 74L, 51L, 62L, 68L, 63L,
        63L, 64L, 53L, 58L, 70L, 62L, 65L, 60L, 65L, 74L, 58L, 54L,
        77L, 47L, 54L, 65L, 57L, 54L, 61L, 65L, 70L, 70L, 56L, 60L,
        62L, 75L, 54L, 62L, 71L, 50L, 65L, 62L, 41L, 54L, 54L, 57L,
        65L, 69L, 46L, 58L, 48L, 46L, 54L, 60L, 60L, 48L, 54L, 63L,
        69L, 56L, 61L, 53L, 63L, 77L, 64L, 58L, 70L, 44L, 56L, 60L,
        64L, 82L, 68L, 52L, 58L, 73L, 49L, 71L, 70L, 53L, 78L, 55L,
        84L, 59L, 72L, 71L, 61L, 62L, 65L, 69L, 48L, 58L, 53L, 67L,
        74L, 55L, 73L, 76L, 57L, 52L, 68L, 58L, 58L, 69L, 50L, 52L,
        62L, 53L, 69L, 56L, 79L, 71L, 56L, 61L, 65L, 55L, 61L, 58L,
        62L, 58L, 66L, 68L, 69L, 66L, 45L, 65L, 51L, 59L, 61L, 70L,
        58L, 60L, 69L, 54L, 79L, 69L, 69L, 58L, 57L, 67L, 52L, 65L,
        71L, 62L, 66L, 53L, 76L, 59L, 68L, 49L, 57L, 80L, 75L, 64L,
        51L, 64L, 71L, 58L, 68L, 58L, 65L, 53L, 69L, 61L, 57L, 68L,
        65L, 69L, 54L, 45L, 76L, 48L, 57L, 70L, 43L, 56L, 56L, 59L,
        57L, 71L, 59L, 64L, 60L, 51L, 58L, 57L, 55L, 72L, 59L, 56L,
        65L, 56L, 46L, 64L, 51L, 59L, 61L, 57L, 60L, 56L, 71L, 62L,
        61L, 60L, 66L, 59L, 57L, 60L, 60L, 66L, 71L, 70L, 60L, 55L,
        62L, 54L, 60L, 68L, 53L, 67L, 75L, 71L, 59L, 77L, 65L, 53L,
        67L, 61L, 73L, 73L, 70L, 59L, 57L, 60L, 66L, 57L, 65L, 60L,
        69L, 56L, 70L, 59L, 69L, 61L, 76L, 61L, 53L, 65L, 60L, 57L,
        57L, 61L, 47L, 71L, 54L, 58L, 58L, 67L, 67L, 58L, 57L, 60L,
        53L, 59L, 66L, 70L, 72L, 54L, 59L, 59L, 52L, 63L, 57L, 60L,
        64L, 56L, 51L, 69L, 65L, 68L, 61L, 51L, 73L, 55L, 56L, 77L,
        63L, 61L, 68L, 62L, 71L, 65L, 67L, 70L, 56L, 64L, 66L, 55L,
        65L, 65L, 59L, 57L, 50L, 64L, 67L, 62L, 67L, 70L, 74L, 68L,
        69L, 64L, 54L, 51L, 52L, 54L, 56L, 71L, 76L, 61L, 48L, 62L,
        66L, 66L, 71L, 60L, 65L, 67L, 66L, 67L, 72L, 61L, 69L, 65L,
        63L, 69L, 75L, 62L, 54L, 71L, 67L, 62L, 52L, 68L, 66L, 67L,
        65L, 71L, 54L, 69L, 69L, 52L, 57L, 59L, 74L, 61L, 57L, 63L,
        66L, 56L, 61L, 51L, 67L, 72L, 63L, 63L, 50L, 66L, 71L, 65L,
        60L, 62L, 58L, 62L, 68L, 56L, 51L, 70L, 51L, 66L, 61L, 54L,
        60L, 56L, 64L, 57L, 63L, 65L, 70L, 57L, 59L, 56L, 62L, 65L,
        59L, 51L, 63L, 65L, 52L, 62L, 55L, 60L, 62L, 64L, 51L, 85L,
        65L, 69L, 57L, 61L, 48L, 81L, 72L, 60L, 59L, 69L, 59L, 61L,
        56L, 65L, 61L, 75L, 64L, 69L, 49L, 60L, 64L, 54L, 71L, 59L,
        62L, 61L, 64L, 55L, 69L, 63L, 55L, 72L, 66L, 68L, 56L, 57L,
        71L, 54L, 70L, 63L, 65L, 59L, 69L, 55L, 71L, 65L, 53L, 69L,
        58L, 74L, 54L, 62L, 51L, 46L, 51L, 53L, 69L, 71L, 68L, 66L,
        56L, 66L, 69L, 59L, 60L, 63L, 58L, 67L, 66L, 69L, 53L, 60L,
        49L, 56L, 55L, 66L, 60L, 51L, 66L, 64L, 43L, 66L, 61L, 68L,
        69L, 71L, 77L, 61L, 45L, 53L, 56L, 54L, 67L, 57L, 65L, 68L,
        75L, 65L, 63L, 67L, 63L, 60L, 63L, 59L, 45L, 65L, 85L, 63L,
        60L, 61L, 69L, 68L, 64L, 73L, 68L, 58L, 56L, 57L, 50L, 55L,
        70L, 49L, 66L, 67L, 76L, 65L, 64L, 68L, 62L, 60L, 61L, 64L,
        66L, 63L, 67L, 68L, 67L, 59L, 49L, 52L, 62L, 59L, 60L, 66L,
        57L, 56L, 74L, 65L, 75L, 69L, 64L, 74L, 59L, 69L, 76L, 65L,
        74L, 66L, 54L, 73L, 61L, 69L, 72L, 52L, 69L, 63L, 61L, 67L,
        60L, 50L, 54L, 80L, 73L, 67L, 53L, 68L, 61L, 60L, 53L, 66L,
        74L, 60L, 60L, 62L, 56L, 59L, 70L, 71L, 84L, 74L, 70L, 68L,
        68L, 61L, 59L, 57L, 53L, 60L, 55L, 57L, 64L, 47L, 77L, 71L,
        70L, 57L, 54L, 62L, 56L, 63L, 70L, 52L, 58L, 67L, 58L, 64L,
        61L, 65L, 56L, 70L, 75L, 68L, 68L, 61L, 48L, 51L, 49L, 61L,
        62L, 67L, 61L, 51L, 55L, 59L, 61L, 68L, 61L, 69L, 72L, 69L,
        58L, 58L, 73L, 58L, 53L, 54L, 58L, 66L, 71L, 72L, 62L, 63L,
        70L, 74L, 72L, 70L, 63L, 75L, 54L, 63L, 54L, 59L, 55L, 55L,
        66L, 60L, 62L, 67L, 55L, 63L, 61L, 69L, 63L, 56L, 58L, 60L,
        56L, 66L, 44L, 58L, 62L, 54L, 71L, 52L, 70L, 76L, 51L, 57L,
        63L, 66L, 58L, 64L, 67L, 55L, 65L, 62L, 51L, 49L, 66L, 72L,
        69L, 64L, 65L, 40L, 68L, 59L, 56L, 70L, 62L, 64L, 68L, 61L,
        65L, 53L, 67L, 71L, 69L, 58L, 57L, 60L, 55L, 74L, 53L, 75L,
        63L, 57L, 75L, 60L, 51L, 59L, 42L, 74L, 90L, 72L, 62L, 62L,
        55L, 64L, 56L, 59L, 61L, 82L, 52L, 51L, 64L, 63L, 65L, 73L,
        49L, 56L, 63L, 69L, 59L, 83L, 72L, 57L, 56L, 62L, 45L, 69L,
        64L, 63L, 61L, 52L, 50L, 71L, 70L, 60L, 54L, 59L, 71L, 64L,
        68L, 66L, 73L, 57L, 64L, 57L, 66L, 59L, 71L, 69L, 65L, 54L,
        61L, 52L, 57L, 67L, 62L, 73L, 53L, 57L, 64L, 60L, 66L, 53L,
        52L, 51L, 63L, 45L, 68L, 67L, 62L, 64L, 57L, 64L, 69L, 58L,
        55L, 46L, 67L, 60L, 70L, 48L, 57L, 48L, 52L, 61L, 63L, 62L,
        63L, 60L, 65L, 69L, 63L, 82L, 65L, 64L, 56L, 64L, 70L, 54L,
        54L, 51L, 58L, 59L, 64L, 56L, 59L, 57L, 65L, 63L, 63L, 77L,
        63L, 57L, 71L, 66L, 52L, 56L, 61L, 74L, 56L, 64L, 73L, 46L,
        57L, 79L, 60L, 62L, 69L, 47L, 69L, 66L, 62L, 46L, 64L, 73L,
        52L, 55L, 53L, 59L, 56L, 67L, 70L, 58L, 60L, 63L, 65L, 48L,
        68L, 59L, 55L, 49L, 62L, 62L, 50L, 55L, 57L, 54L, 57L, 62L,
        53L, 57L, 59L, 59L, 56L, 62L, 63L, 57L, 61L, 65L, 74L, 52L,
        62L, 64L, 54L, 57L, 62L, 61L, 63L, 50L, 65L, 51L, 64L, 43L,
        65L, 63L, 66L, 68L, 72L, 63L, 61L, 59L, 64L, 60L, 60L, 54L,
        65L, 62L, 60L, 68L, 64L, 55L, 60L, 72L, 66L, 65L, 52L, 57L,
        73L, 55L, 54L, 56L, 61L, 66L, 57L, 60L, 62L, 61L, 76L, 77L,
        66L, 70L, 71L, 79L, 72L, 55L, 66L, 48L, 52L, 61L, 64L, 66L,
        65L, 54L, 53L, 64L, 46L, 68L, 57L, 70L, 60L, 63L, 67L, 64L,
        71L, 64L, 68L, 66L, 59L, 57L, 54L, 60L, 73L, 56L, 50L, 47L,
        69L, 67L, 51L, 75L, 61L, 65L, 63L, 54L, 66L, 69L, 82L, 70L,
        60L, 59L, 66L, 57L, 68L, 62L, 70L, 52L, 66L, 67L, 64L, 46L,
        63L, 61L, 60L, 55L, 76L, 60L, 62L, 61L, 69L, 65L, 78L, 58L,
        59L, 56L, 66L, 60L, 56L, 71L, 65L, 75L, 65L, 54L, 50L, 59L,
        69L, 63L, 64L, 57L, 56L, 59L, 66L, 65L, 58L, 70L, 53L, 61L,
        46L, 61L, 65L, 69L, 67L, 63L, 54L, 52L, 66L, 57L, 63L, 51L,
        56L, 50L, 57L, 59L, 54L, 68L, 53L, 67L, 72L, 70L, 63L, 52L,
        47L, 62L, 62L, 66L, 56L, 51L, 71L, 67L, 61L, 44L, 68L, 63L,
        62L, 75L, 65L, 76L, 49L, 70L, 59L, 72L, 73L, 60L, 64L, 48L,
        70L, 60L, 53L, 67L, 71L, 71L, 52L, 69L, 79L, 59L, 53L, 74L,
        69L, 62L, 55L, 66L, 63L, 50L, 73L, 51L, 68L, 68L, 57L, 61L,
        63L, 59L, 74L, 61L, 58L, 76L, 52L, 59L, 55L, 57L, 61L, 64L,
        77L, 65L, 57L, 55L, 60L, 68L, 76L, 67L, 64L, 50L, 70L, 61L,
        61L, 71L, 55L, 70L, 64L, 68L, 64L, 52L, 65L, 67L, 64L, 56L,
        59L, 63L, 78L, 58L, 55L, 61L, 61L, 73L, 53L, 56L, 49L, 55L,
        52L, 62L, 75L, 62L, 68L, 57L, 48L, 56L, 43L, 60L, 60L, 71L,
        79L, 67L, 78L, 74L, 54L, 63L, 67L, 69L, 67L, 67L, 58L, 68L,
        57L, 61L, 55L, 60L, 75L, 72L, 79L, 67L, 63L, 63L, 53L, 57L,
        65L, 62L, 58L, 59L, 73L, 49L, 55L, 70L, 53L, 70L, 55L, 55L,
        71L, 76L, 72L, 60L, 47L, 71L, 67L, 75L, 57L, 68L, 70L, 75L,
        65L, 67L, 77L, 80L, 80L, 53L, 57L, 66L, 66L, 62L, 51L, 58L,
        56L, 54L, 61L, 68L, 57L, 56L, 50L, 72L, 53L, 58L, 53L, 49L,
        64L, 66L, 43L, 75L, 67L, 63L, 45L, 63L, 46L, 72L, 61L, 61L,
        70L, 51L, 76L, 61L)
stopifnot(identical(indegree(model), id))

od <- c(61L, 60L, 64L, 64L, 67L, 66L, 59L, 58L, 58L, 62L, 60L, 55L,
        64L, 47L, 65L, 50L, 53L, 70L, 78L, 67L, 60L, 62L, 65L, 62L,
        69L, 58L, 56L, 60L, 74L, 63L, 60L, 54L, 69L, 80L, 61L, 53L,
        60L, 70L, 65L, 56L, 64L, 73L, 69L, 62L, 59L, 70L, 65L, 59L,
        54L, 56L, 68L, 71L, 65L, 76L, 63L, 81L, 60L, 57L, 53L, 58L,
        59L, 67L, 71L, 74L, 66L, 68L, 59L, 69L, 58L, 54L, 58L, 71L,
        65L, 55L, 48L, 62L, 62L, 58L, 78L, 51L, 67L, 63L, 68L, 60L,
        65L, 56L, 67L, 68L, 70L, 47L, 59L, 69L, 66L, 75L, 62L, 58L,
        58L, 80L, 60L, 71L, 69L, 74L, 59L, 58L, 59L, 61L, 70L, 62L,
        66L, 51L, 54L, 61L, 65L, 61L, 63L, 59L, 74L, 50L, 50L, 74L,
        69L, 67L, 56L, 69L, 58L, 55L, 64L, 58L, 47L, 66L, 56L, 52L,
        59L, 62L, 64L, 58L, 61L, 78L, 75L, 81L, 61L, 69L, 64L, 74L,
        59L, 62L, 64L, 58L, 65L, 59L, 62L, 49L, 48L, 46L, 67L, 65L,
        75L, 69L, 71L, 60L, 73L, 66L, 51L, 67L, 76L, 57L, 72L, 82L,
        58L, 73L, 54L, 52L, 71L, 66L, 68L, 53L, 75L, 63L, 67L, 63L,
        55L, 63L, 39L, 61L, 72L, 66L, 61L, 72L, 56L, 61L, 59L, 51L,
        80L, 56L, 60L, 64L, 55L, 68L, 69L, 66L, 58L, 65L, 68L, 54L,
        61L, 66L, 65L, 64L, 63L, 69L, 64L, 62L, 65L, 58L, 62L, 58L,
        56L, 57L, 68L, 65L, 53L, 60L, 72L, 62L, 74L, 61L, 66L, 62L,
        59L, 61L, 63L, 66L, 65L, 57L, 64L, 53L, 57L, 57L, 56L, 59L,
        52L, 48L, 59L, 62L, 72L, 65L, 54L, 68L, 59L, 45L, 52L, 53L,
        58L, 66L, 52L, 64L, 68L, 57L, 63L, 65L, 67L, 76L, 64L, 59L,
        73L, 72L, 58L, 55L, 72L, 54L, 65L, 56L, 63L, 61L, 59L, 63L,
        65L, 63L, 85L, 63L, 59L, 66L, 51L, 66L, 61L, 50L, 78L, 43L,
        63L, 63L, 63L, 56L, 63L, 71L, 78L, 56L, 59L, 72L, 43L, 57L,
        74L, 67L, 58L, 48L, 49L, 61L, 51L, 54L, 69L, 61L, 61L, 69L,
        57L, 62L, 53L, 66L, 53L, 64L, 70L, 67L, 66L, 72L, 66L, 67L,
        68L, 61L, 63L, 72L, 59L, 68L, 55L, 53L, 67L, 63L, 68L, 71L,
        48L, 70L, 44L, 67L, 62L, 48L, 60L, 75L, 60L, 54L, 57L, 73L,
        74L, 67L, 54L, 53L, 68L, 71L, 63L, 75L, 68L, 55L, 66L, 67L,
        58L, 53L, 57L, 66L, 61L, 66L, 55L, 50L, 65L, 59L, 63L, 67L,
        65L, 70L, 57L, 49L, 64L, 59L, 62L, 68L, 74L, 56L, 66L, 63L,
        66L, 71L, 70L, 64L, 61L, 49L, 50L, 71L, 70L, 60L, 70L, 58L,
        66L, 55L, 62L, 63L, 60L, 55L, 51L, 57L, 48L, 48L, 60L, 65L,
        69L, 52L, 80L, 73L, 69L, 72L, 63L, 65L, 59L, 56L, 51L, 69L,
        55L, 69L, 54L, 66L, 69L, 77L, 68L, 63L, 66L, 61L, 64L, 70L,
        48L, 52L, 59L, 62L, 55L, 58L, 68L, 59L, 79L, 81L, 60L, 66L,
        53L, 73L, 54L, 63L, 85L, 68L, 65L, 77L, 65L, 50L, 58L, 55L,
        57L, 66L, 81L, 57L, 63L, 74L, 60L, 59L, 58L, 60L, 57L, 58L,
        56L, 61L, 64L, 60L, 68L, 68L, 62L, 57L, 48L, 51L, 51L, 58L,
        73L, 53L, 50L, 62L, 57L, 55L, 56L, 62L, 55L, 60L, 70L, 62L,
        59L, 60L, 70L, 57L, 64L, 53L, 65L, 65L, 61L, 59L, 65L, 60L,
        60L, 56L, 57L, 52L, 58L, 66L, 51L, 55L, 56L, 60L, 54L, 59L,
        68L, 54L, 72L, 70L, 57L, 65L, 66L, 63L, 66L, 71L, 57L, 69L,
        69L, 57L, 68L, 57L, 67L, 71L, 50L, 63L, 75L, 64L, 63L, 51L,
        57L, 60L, 63L, 65L, 76L, 76L, 55L, 59L, 66L, 56L, 59L, 68L,
        79L, 61L, 79L, 49L, 66L, 64L, 62L, 59L, 60L, 68L, 65L, 81L,
        56L, 56L, 72L, 59L, 73L, 43L, 66L, 60L, 67L, 60L, 65L, 72L,
        60L, 73L, 72L, 56L, 56L, 71L, 62L, 78L, 70L, 53L, 61L, 77L,
        66L, 71L, 59L, 50L, 65L, 63L, 70L, 66L, 67L, 48L, 45L, 66L,
        63L, 68L, 70L, 63L, 56L, 48L, 67L, 62L, 67L, 53L, 69L, 60L,
        64L, 70L, 51L, 63L, 72L, 51L, 59L, 74L, 57L, 63L, 61L, 60L,
        57L, 57L, 68L, 54L, 67L, 59L, 47L, 63L, 53L, 69L, 55L, 71L,
        56L, 73L, 57L, 70L, 70L, 41L, 66L, 64L, 55L, 60L, 78L, 63L,
        66L, 72L, 70L, 67L, 70L, 61L, 75L, 53L, 59L, 66L, 70L, 67L,
        57L, 60L, 51L, 54L, 58L, 56L, 61L, 57L, 80L, 68L, 70L, 53L,
        62L, 59L, 68L, 69L, 53L, 68L, 67L, 61L, 52L, 61L, 61L, 64L,
        63L, 63L, 60L, 55L, 56L, 63L, 61L, 72L, 54L, 54L, 47L, 68L,
        64L, 67L, 68L, 69L, 52L, 59L, 56L, 61L, 66L, 50L, 66L, 67L,
        57L, 55L, 54L, 48L, 62L, 75L, 70L, 61L, 55L, 63L, 59L, 65L,
        75L, 65L, 65L, 71L, 69L, 48L, 57L, 64L, 61L, 58L, 68L, 69L,
        78L, 63L, 68L, 59L, 66L, 68L, 60L, 55L, 66L, 59L, 75L, 66L,
        48L, 67L, 56L, 60L, 62L, 48L, 66L, 54L, 72L, 53L, 55L, 61L,
        67L, 59L, 68L, 70L, 53L, 53L, 78L, 50L, 62L, 70L, 62L, 67L,
        67L, 72L, 56L, 55L, 62L, 71L, 48L, 72L, 56L, 55L, 58L, 63L,
        63L, 73L, 65L, 57L, 60L, 70L, 69L, 50L, 58L, 67L, 60L, 68L,
        63L, 48L, 63L, 77L, 61L, 58L, 87L, 66L, 61L, 73L, 69L, 65L,
        60L, 69L, 44L, 61L, 46L, 46L, 63L, 59L, 47L, 58L, 57L, 61L,
        79L, 56L, 60L, 57L, 64L, 55L, 71L, 62L, 65L, 62L, 69L, 64L,
        56L, 61L, 56L, 55L, 53L, 56L, 59L, 72L, 68L, 47L, 74L, 62L,
        53L, 54L, 56L, 66L, 65L, 65L, 69L, 61L, 72L, 58L, 44L, 58L,
        56L, 55L, 67L, 68L, 56L, 65L, 65L, 62L, 66L, 67L, 57L, 89L,
        69L, 77L, 79L, 77L, 53L, 66L, 49L, 63L, 60L, 62L, 69L, 61L,
        68L, 70L, 55L, 50L, 62L, 68L, 62L, 63L, 51L, 51L, 78L, 66L,
        58L, 76L, 68L, 67L, 63L, 55L, 63L, 66L, 64L, 62L, 59L, 57L,
        59L, 59L, 66L, 61L, 58L, 57L, 71L, 56L, 69L, 64L, 55L, 59L,
        76L, 71L, 74L, 73L, 60L, 46L, 72L, 65L, 59L, 52L, 76L, 54L,
        59L, 70L, 64L, 59L, 57L, 66L, 68L, 62L, 57L, 56L, 51L, 63L,
        71L, 61L, 69L, 60L, 66L, 61L, 54L, 48L, 59L, 61L, 49L, 50L,
        66L, 77L, 63L, 57L, 48L, 51L, 71L, 57L, 61L, 65L, 58L, 60L,
        70L, 67L, 56L, 63L, 56L, 47L, 62L, 75L, 68L, 72L, 82L, 67L,
        70L, 36L, 70L, 63L, 63L, 63L, 58L, 66L, 60L, 59L, 65L, 68L,
        61L, 70L, 65L, 50L, 65L, 52L, 65L, 47L, 68L, 68L, 60L, 56L,
        52L, 60L, 61L, 59L, 66L, 51L, 68L, 72L, 85L, 64L, 49L, 48L,
        54L, 45L, 76L, 67L, 56L, 57L, 62L, 68L, 67L, 74L, 67L, 53L,
        63L, 64L, 54L, 61L, 62L, 61L, 53L, 55L, 61L, 53L, 51L, 59L,
        50L, 68L, 53L, 55L, 78L, 63L, 62L, 62L, 58L, 59L, 58L, 54L,
        68L, 61L, 74L, 76L, 67L, 67L, 66L, 48L, 77L, 64L, 56L, 75L,
        55L, 85L, 65L, 64L, 59L, 65L, 65L, 55L, 53L, 61L, 58L, 53L,
        59L, 45L, 67L, 58L, 50L, 63L, 52L, 61L, 75L, 70L, 57L, 66L,
        67L, 73L, 79L, 67L, 69L, 54L, 51L, 73L, 61L, 78L, 55L, 74L,
        55L, 57L, 64L, 52L, 61L, 59L, 53L, 73L, 57L, 53L, 55L, 59L,
        65L, 60L, 58L, 55L, 70L, 49L, 72L, 60L, 49L, 61L, 52L, 71L,
        59L, 61L, 68L, 66L, 53L, 63L, 69L, 55L, 62L, 52L, 84L, 59L,
        62L, 55L, 73L, 56L, 60L, 57L, 64L, 56L, 68L, 53L, 68L, 49L,
        70L, 51L, 54L, 62L, 57L, 63L, 69L, 53L, 60L, 63L, 53L, 80L,
        60L, 68L, 57L, 59L, 59L, 64L, 53L, 62L, 68L, 62L, 58L, 55L,
        67L, 50L, 58L, 54L, 57L, 68L, 68L, 60L, 60L, 74L, 56L, 75L,
        60L, 71L, 55L, 68L, 62L, 59L, 60L, 65L, 57L, 54L, 61L, 66L,
        55L, 48L, 44L, 61L, 67L, 72L, 66L, 66L, 59L, 62L, 51L, 58L,
        55L, 68L, 68L, 61L, 65L, 64L, 55L, 80L, 59L, 58L, 64L, 74L,
        72L, 64L, 68L, 70L, 67L, 74L, 63L, 65L, 56L, 65L, 49L, 53L,
        65L, 64L, 64L, 79L, 72L, 66L, 56L, 61L, 48L, 63L, 57L, 77L,
        54L, 69L, 75L, 58L, 56L, 71L, 59L, 65L, 58L, 65L, 60L, 49L,
        56L, 64L, 51L, 69L, 69L, 59L, 67L, 74L, 60L, 73L, 60L, 58L,
        50L, 56L, 70L, 57L, 65L, 67L, 57L, 58L, 54L, 65L, 66L, 61L,
        68L, 58L, 74L, 65L, 62L, 66L, 61L, 66L, 70L, 66L, 76L, 59L,
        61L, 60L, 64L, 65L, 55L, 63L, 70L, 62L, 58L, 55L, 72L, 61L,
        68L, 47L, 50L, 69L, 59L, 67L, 76L, 74L, 68L, 55L, 64L, 55L,
        64L, 72L, 56L, 69L, 66L, 60L, 52L, 54L, 60L, 61L, 61L, 59L,
        62L, 73L, 62L, 61L, 73L, 59L, 54L, 51L, 67L, 63L, 57L, 62L,
        58L, 68L, 69L, 63L, 73L, 55L, 83L, 59L, 44L, 67L, 62L, 60L,
        69L, 52L, 63L, 63L, 58L, 52L, 67L, 63L, 61L, 52L, 67L, 63L,
        59L, 60L, 70L, 64L, 50L, 53L, 77L, 73L, 68L, 60L, 46L, 67L,
        71L, 68L, 75L, 63L, 62L, 74L, 54L, 67L, 62L, 60L, 69L, 56L,
        67L, 51L, 67L, 59L, 62L, 61L, 71L, 58L, 58L, 60L, 62L, 79L,
        64L, 52L, 68L, 56L, 56L, 59L, 55L, 58L, 54L, 67L, 54L, 64L,
        56L, 80L, 72L, 72L, 63L, 58L, 60L, 59L, 58L, 57L, 54L, 68L,
        58L, 64L, 50L, 60L, 56L, 60L, 67L, 63L, 50L, 59L, 52L, 70L,
        55L, 60L, 81L, 63L, 62L, 56L, 66L, 54L, 53L, 69L, 55L, 56L,
        66L, 58L, 70L, 60L, 76L, 58L, 58L, 78L, 60L, 67L, 58L, 54L,
        62L, 68L, 57L, 44L, 54L, 50L, 70L, 54L, 73L, 56L, 73L, 72L,
        68L, 73L, 68L, 55L, 76L, 56L, 58L, 66L, 51L, 56L, 66L, 62L,
        59L, 60L, 65L, 55L, 59L, 69L, 63L, 55L, 57L, 70L, 56L, 73L,
        58L, 78L, 72L, 57L, 70L, 62L, 70L, 51L, 58L, 80L, 51L, 56L,
        56L, 62L, 73L, 59L, 64L, 51L, 62L, 66L, 54L, 68L, 59L, 47L,
        57L, 68L, 69L, 58L, 61L, 62L, 54L, 66L, 63L, 68L, 70L, 53L,
        59L, 59L, 47L, 63L, 68L, 65L, 66L, 55L, 61L, 65L, 65L, 64L,
        51L, 59L, 72L, 58L, 67L, 67L, 64L, 68L, 43L, 57L, 47L, 56L,
        47L, 71L, 57L, 50L, 74L, 61L, 55L, 67L, 68L, 50L, 75L, 64L,
        62L, 61L, 60L, 63L, 65L, 58L, 64L, 74L, 69L, 82L, 52L, 64L,
        71L, 54L, 56L, 67L, 69L, 59L, 61L, 55L, 49L, 62L, 49L, 68L,
        72L, 82L, 48L, 68L, 50L, 54L, 70L, 74L, 53L, 79L, 56L, 55L,
        52L, 67L, 51L, 64L)
stopifnot(identical(outdegree(model), od))

## Check that initialisation raises an invalid rate error.
model <- SIR(u0 = data.frame(S = rep(99, 2), I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = -0.077)

res <- assertError(run(model, solver = "ssm"))
check_error(res, "Invalid rate detected (non-finite or < 0.0).")

res <- assertError(run(model, solver = "aem"))
check_error(res, "Invalid rate detected (non-finite or < 0.0).")

## Check that an invalid solver argument raises an error.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

res <- assertError(.Call(SimInf:::SIR_run, model, c("ssm", "ssm")))
check_error(res, "Invalid 'solver' value.")
res <- assertError(.Call(SimInf:::SIR_run, model, "non-existing-solver"))
check_error(res, "Invalid 'solver' value.")
res <- assertError(.Call(SimInf:::SIR_run, model, NA_character_))
check_error(res, "Invalid 'solver' value.")
res <- assertError(.Call(SimInf:::SIR_run, model, 5))
check_error(res, "Invalid 'solver' value.")

## Trigger a negative state error
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0)
model@S@x[1] <- -1000
set.seed(123)
res <- assertError(.Call(SimInf:::SIR_run, model, "ssm"))
check_error(res, "Negative state detected.")
res <- assertError(.Call(SimInf:::SIR_run, model, "aem"))
check_error(res, "Negative state detected.")
