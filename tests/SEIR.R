## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2020 Stefan Widgren
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
res <- assertError(SEIR(u0 = "u0"))
check_error(res, "Missing columns in u0.")

u0 <- data.frame(S  = c(0, 1, 2, 3, 4, 5),
                 E  = c(0, 0, 0, 0, 0, 0),
                 I  = c(0, 0, 0, 0, 0, 0),
                 R  = c(0, 0, 0, 0, 0, 0))

## Check missing columns in u0
res <- assertError(SEIR(u0 = u0[, c("E", "I", "R"), drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- assertError(SEIR(u0 = u0[, c("S", "I", "R"), drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- assertError(SEIR(u0 = u0[, c("S", "E", "R"), drop = FALSE]))
check_error(res, "Missing columns in u0.")

res <- assertError(SEIR(u0 = u0[, c("S", "E", "I"), drop = FALSE]))
check_error(res, "Missing columns in u0.")

## Check missing beta
res <- assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               epsilon = 0.5,
                               gamma   = 0.5))
check_error(res, "'beta' must be numeric of length 1 or 'nrow(u0)'.")

## Check missing epsilon
res <- assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               gamma   = 0.5))
check_error(res, "'epsilon' must be numeric of length 1 or 'nrow(u0)'.")

## Check missing gamma
res <- assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = 0.5))
check_error(res, "'gamma' must be numeric of length 1 or 'nrow(u0)'.")

## Check non-numeric beta
res <- assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = "0.5",
                               epsilon = 0.3,
                               gamma   = 0.1))
check_error(res, "'beta' must be numeric of length 1 or 'nrow(u0)'.")

## Check non-numeric epsilon
res <- assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = "0.3",
                               gamma   = 0.1))
check_error(res, "'epsilon' must be numeric of length 1 or 'nrow(u0)'.")

## Check non-numeric gamma
res <- assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = 0.3,
                               gamma   = "0.1"))
check_error(res, "'gamma' must be numeric of length 1 or 'nrow(u0)'.")

## Check that length of beta equals 1
res <- assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = c(0.5, 0.5),
                               epsilon = 0.3,
                               gamma   = 0.1))
check_error(res, "'beta' must be numeric of length 1 or 'nrow(u0)'.")

## Check that length of epsilon equals 1
res <- assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = c(0.3, 0.3),
                               gamma   = 0.1))
check_error(res, "'epsilon' must be numeric of length 1 or 'nrow(u0)'.")

## Check that length of gamma equals 1
res <- assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = 0.3,
                               gamma   = c(0.1, 0.1)))
check_error(res, "'gamma' must be numeric of length 1 or 'nrow(u0)'.")

## Check extraction of data from 'suscpetible', 'infected' and
## 'recovered' compartments
model <- SEIR(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              beta    = 0,
              epsilon = 0,
              gamma   = 0)

result <- run(model)

S_expected <- structure(c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L,
                          1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
                          2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L,
                          3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
                          4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                        .Dim = c(6L, 10L))

S_observed <- trajectory(result, compartments = "S", as.is = TRUE)
stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L))

I_observed <- trajectory(result, compartments = "I", as.is = TRUE)
stopifnot(identical(I_observed, I_expected))

R_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L))

R_observed <- trajectory(result, compartments = "R", as.is = TRUE)
stopifnot(identical(R_observed, R_expected))

R_expected <- data.frame(
    node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L,
             3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L,
             1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L),
    time = c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L,
             1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L,
             4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 7L,
             7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L,
             9L),
    R = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L))
R_observed <- trajectory(result, compartments = "R")
stopifnot(identical(R_observed, R_expected))

## Extract the number of recovered individuals in the first node after
## each time step in the simulation
R_expected <- data.frame(
    node = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    time = 0:9L,
    R = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L))
R_observed <- trajectory(result, compartments = "R", index = 1)
stopifnot(identical(R_observed, R_expected))

R_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                       .Dim = c(1L, 10L))
R_observed <- trajectory(result, compartments = "R", index = 1, as.is = TRUE)
stopifnot(identical(R_observed, R_expected))

## Extract the number of recovered individuals in the first and third
## node after each time step in the simulation
R_expected <- data.frame(
    node = c(1L, 3L, 1L, 3L, 1L, 3L, 1L, 3L, 1L, 3L,
             1L, 3L, 1L, 3L, 1L, 3L, 1L, 3L, 1L, 3L),
    time = c(0L, 0L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L,
             5L, 5L, 6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L),
    R = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L))
R_observed <- trajectory(result, compartments = "R", index = c(1, 3))
stopifnot(identical(R_observed, R_expected))

R_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(2L, 10L))
R_observed <- trajectory(result, compartments = "R", index = c(1, 3),
                         as.is = TRUE)
stopifnot(identical(R_observed, R_expected))

## A more complex test to extract data from U from a trajectory of 6
## nodes and 10 time steps (0:9). Generate events that adds one
## individual per node in each compartment each day. To achieve this,
## using the template SEIR model, first add them to the 'S'
## compartment, and then modify the select matrix E and update the
## select index.
events <- do.call("rbind", lapply(seq_len(9), function(time) {
    do.call("rbind", lapply(seq_len(6), function(node) {
        data.frame(event = 1, time  = time, node  = node, dest  = 0,
                   n = 1, proportion = 0, select = rep(1, 4), shift = 0)
    }))
}))
model <- SEIR(u0 = data.frame(S = c(110, 210, 310, 410, 510, 610),
                              E = c(120, 220, 320, 420, 520, 620),
                              I = c(130, 230, 330, 430, 530, 630),
                              R = c(140, 240, 340, 440, 540, 640)),
              tspan   = 0:9, events  = events, beta    = 0,
              epsilon = 0, gamma   = 0)
model@events@E <- as(diag(4), "dgCMatrix")
model@events@select <- rep(1:4, length.out = length(model@events@select))

# Check that this fails because rownames (compartments) are missing
res <- assertError(run(model))
check_error(
    res,
    "'S' and 'E' must have rownames matching the compartments.",
    FALSE)

rownames(model@events@E) <- c("S", "E", "I", "R")
result <- run(model)

U_expected <- data.frame(
    node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L,
             3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L,
             1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L),
    time = c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L,
             1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L,
             4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 7L,
             7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L,
             9L),
    S = c(110L, 210L, 310L, 410L, 510L, 610L, 111L, 211L, 311L,
          411L, 511L, 611L, 112L, 212L, 312L, 412L, 512L, 612L, 113L, 213L,
          313L, 413L, 513L, 613L, 114L, 214L, 314L, 414L, 514L, 614L, 115L,
          215L, 315L, 415L, 515L, 615L, 116L, 216L, 316L, 416L, 516L, 616L,
          117L, 217L, 317L, 417L, 517L, 617L, 118L, 218L, 318L, 418L, 518L,
          618L, 119L, 219L, 319L, 419L, 519L, 619L),
    E = c(120L, 220L,
          320L, 420L, 520L, 620L, 121L, 221L, 321L, 421L, 521L, 621L, 122L,
          222L, 322L, 422L, 522L, 622L, 123L, 223L, 323L, 423L, 523L, 623L,
          124L, 224L, 324L, 424L, 524L, 624L, 125L, 225L, 325L, 425L, 525L,
          625L, 126L, 226L, 326L, 426L, 526L, 626L, 127L, 227L, 327L, 427L,
          527L, 627L, 128L, 228L, 328L, 428L, 528L, 628L, 129L, 229L, 329L,
          429L, 529L, 629L),
    I = c(130L, 230L, 330L, 430L, 530L, 630L,
          131L, 231L, 331L, 431L, 531L, 631L, 132L, 232L, 332L, 432L, 532L,
          632L, 133L, 233L, 333L, 433L, 533L, 633L, 134L, 234L, 334L, 434L,
          534L, 634L, 135L, 235L, 335L, 435L, 535L, 635L, 136L, 236L, 336L,
          436L, 536L, 636L, 137L, 237L, 337L, 437L, 537L, 637L, 138L, 238L,
          338L, 438L, 538L, 638L, 139L, 239L, 339L, 439L, 539L, 639L),
    R = c(140L, 240L, 340L, 440L, 540L, 640L, 141L, 241L, 341L,
          441L, 541L, 641L, 142L, 242L, 342L, 442L, 542L, 642L, 143L,
          243L, 343L, 443L, 543L, 643L, 144L, 244L, 344L, 444L, 544L,
          644L, 145L, 245L, 345L, 445L, 545L, 645L, 146L, 246L, 346L,
          446L, 546L, 646L, 147L, 247L, 347L, 447L, 547L, 647L, 148L,
          248L, 348L, 448L, 548L, 648L, 149L, 249L, 349L, 449L, 549L,
          649L))

U_observed <- trajectory(result)
stopifnot(identical(U_observed, U_expected))
U_observed <- trajectory(result, compartments = c("S", "E", "I", "R"),
                         index = 1:6)
stopifnot(identical(U_observed, U_expected))

U_expected <- data.frame(
    node = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    time = 0:9, S = 110:119, I = 130:139)
U_observed <- trajectory(result, compartments = c("S", "I"), index = 1)
stopifnot(identical(U_observed, U_expected))

U_expected <- data.frame(
    node = c(2L, 4L, 2L, 4L, 2L, 4L, 2L, 4L, 2L, 4L,
             2L, 4L, 2L, 4L, 2L, 4L, 2L, 4L, 2L, 4L),
    time = c(0L, 0L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L,
             5L, 5L, 6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L),
    E = c(220L, 420L, 221L, 421L, 222L, 422L, 223L, 423L, 224L, 424L,
          225L, 425L, 226L, 426L, 227L, 427L, 228L, 428L, 229L, 429L),
    R = c(240L, 440L, 241L, 441L, 242L, 442L, 243L, 443L, 244L, 444L,
          245L, 445L, 246L, 446L, 247L, 447L, 248L, 448L, 249L, 449L))
U_observed <- trajectory(result, compartments = c("E", "R"), index = c(2, 4))
stopifnot(identical(U_observed, U_expected))

U_expected <- structure(
    c(220L, 240L, 420L, 440L, 221L, 241L, 421L, 441L, 222L,
      242L, 422L, 442L, 223L, 243L, 423L, 443L, 224L, 244L, 424L, 444L,
      225L, 245L, 425L, 445L, 226L, 246L, 426L, 446L, 227L, 247L, 427L,
      447L, 228L, 248L, 428L, 448L, 229L, 249L, 429L, 449L),
    .Dim = c(4L, 10L))
U_observed <- trajectory(result, compartments = c("E", "R"),
                         index = c(2, 4), as.is = TRUE)
stopifnot(identical(U_observed, U_expected))

## Check prevalence

## Define a tolerance
tol <- 1e-8

p_expected <- data.frame(
    time = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
    prevalence = c(0.253333333333333, 0.253324468085106, 0.253315649867374,
                   0.253306878306878, 0.253298153034301, 0.253289473684211,
                   0.253280839895013, 0.253272251308901, 0.253263707571802,
                   0.253255208333333))
p_observed <- prevalence(result, I ~ S + E + I + R)
stopifnot(identical(p_observed$time, p_expected$time))
stopifnot(all(abs(p_observed$prevalence - p_expected$prevalence) < tol))
p_observed <- prevalence(result, I~.)
stopifnot(identical(p_observed$time, p_expected$time))
stopifnot(all(abs(p_observed$prevalence - p_expected$prevalence) < tol))

p_expected <- structure(
    list(time = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
         prevalence = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
    .Names = c("time", "prevalence"),
    row.names = c(NA, -10L),
    class = "data.frame")
p_observed <- prevalence(result, I ~ S + E + I + R, type = "nop")
stopifnot(identical(p_observed, p_expected))
p_observed <- prevalence(result, I ~ ., type = "nop")
stopifnot(identical(p_observed, p_expected))

p_expected <- data.frame(
    node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L,
             3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L,
             1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L),
    time = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3,
             3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6,
             7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9),
    prevalence = c(0.26, 0.255555555555556, 0.253846153846154,
                   0.252941176470588, 0.252380952380952, 0.252,
                   0.259920634920635, 0.255530973451327, 0.253834355828221,
                   0.252934272300469, 0.252376425855513, 0.251996805111821,
                   0.259842519685039, 0.255506607929515, 0.253822629969419,
                   0.252927400468384, 0.252371916508539, 0.251993620414673,
                   0.259765625, 0.255482456140351, 0.253810975609756,
                   0.252920560747664, 0.252367424242424, 0.251990445859873,
                   0.25968992248062, 0.255458515283843, 0.253799392097264,
                   0.252913752913753, 0.252362948960302, 0.251987281399046,
                   0.259615384615385, 0.255434782608696, 0.253787878787879,
                   0.252906976744186, 0.252358490566038, 0.251984126984127,
                   0.259541984732824, 0.255411255411255, 0.253776435045317,
                   0.252900232018561, 0.252354048964218, 0.251980982567353,
                   0.259469696969697, 0.255387931034483, 0.253765060240964,
                   0.252893518518519, 0.25234962406015, 0.251977848101266,
                   0.259398496240602, 0.255364806866953, 0.253753753753754,
                   0.252886836027714, 0.25234521575985, 0.251974723538705,
                   0.259328358208955, 0.25534188034188, 0.25374251497006,
                   0.252880184331797, 0.252340823970037, 0.251971608832808))
p_observed <- prevalence(result, I ~ ., type = "wnp")
stopifnot(identical(p_observed$node, p_expected$node))
stopifnot(identical(p_observed$time, p_expected$time))
stopifnot(all(abs(p_observed$prevalence - p_expected$prevalence) < tol))
p_observed <- prevalence(result, I ~ S + E + I + R, type = "wnp")
stopifnot(identical(p_observed$node, p_expected$node))
stopifnot(identical(p_observed$time, p_expected$time))
stopifnot(all(abs(p_observed$prevalence - p_expected$prevalence) < tol))

p_expected <- data.frame(
    node = c(2L, 3L, 2L, 3L, 2L, 3L, 2L, 3L, 2L, 3L, 2L, 3L, 2L, 3L, 2L, 3L,
             2L, 3L, 2L, 3L),
    time = c(0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9),
    prevalence = c(0.255555555555556, 0.253846153846154, 0.255530973451327,
                   0.253834355828221, 0.255506607929515, 0.253822629969419,
                   0.255482456140351, 0.253810975609756, 0.255458515283843,
                   0.253799392097264, 0.255434782608696, 0.253787878787879,
                   0.255411255411255, 0.253776435045317, 0.255387931034483,
                   0.253765060240964, 0.255364806866953, 0.253753753753754,
                   0.25534188034188, 0.25374251497006))
p_observed <- prevalence(result, I~., type = "wnp", i = 2:3)
stopifnot(identical(p_observed$node, p_expected$node))
stopifnot(identical(p_observed$time, p_expected$time))
stopifnot(all(abs(p_observed$prevalence - p_expected$prevalence) < tol))
p_observed <- prevalence(result, I ~ S + E + I + R, type = "wnp", i = 2:3)
stopifnot(identical(p_observed$node, p_expected$node))
stopifnot(identical(p_observed$time, p_expected$time))
stopifnot(all(abs(p_observed$prevalence - p_expected$prevalence) < tol))

## Check 'V'
res <- assertError(trajectory(result, "phi"))
check_error(res, "Non-existing compartment(s) in model: 'phi'.")

## Check data
stopifnot(identical(nrow(events_SEIR()), 466692L))
stopifnot(identical(nrow(u0_SEIR()), 1600L))

## Try to plot non-extisting compartment.
res <- assertError(plot(result, compartments = "X"))
check_error(res, "'compartments' must exist in the model.")

## Try to plot with invalid range argument.
res <- assertError(plot(result, range = 1.2))
check_error(res, "'range' must be FALSE or a value between 0 and 1.")

## Check SEIR plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SEIR plot method with range = FALSE
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result, compartments = "S", lty = 1, range = FALSE)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SEIR boxplot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
boxplot(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SEIR pairs plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
pairs(result)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SEIR events plot with no events
model <- SEIR(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              beta    = 0,
              epsilon = 0,
              gamma   = 0)
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(model@events)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check SEIR events plot method
model <- SEIR(u0      = u0_SEIR(),
              tspan   = seq_len(365 * 4),
              events  = events_SEIR(),
              beta    = 0,
              epsilon = 0,
              gamma   = 0)
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(model@events)
dev.off()
stopifnot(file.exists(pdf_file))
unlink(pdf_file)

## Check that C SEIR run function fails for a misspecified SEIR model
res <- assertError(.Call(SimInf:::SEIR_run, NULL, NULL))
check_error(res, "Invalid model.")

res <- assertError(.Call(SimInf:::SEIR_run, "SEIR", NULL))
check_error(res, "Invalid model.")

## Check that an invalid rate error is raised during the simulation.
model <- SEIR(u0 = data.frame(S = rep(100, 2), E = 0, I = 10, R = 0),
              tspan = 1:100,
              beta = 0.16,
              epsilon = -0.3,
              gamma = 0.077)
set.seed(1)
res <- assertError(run(model, solver = "ssm"))
check_error(res, "Invalid rate detected (non-finite or < 0.0).")
res <- assertError(run(model, solver = "aem"))
check_error(res, "Invalid rate detected (non-finite or < 0.0).")
