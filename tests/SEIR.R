## SimInf, a framework for stochastic disease spread simulations
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

## For debugging
sessionInfo()

## Check invalid u0
res <- tools::assertError(SEIR(u0 = "u0"))
stopifnot(length(grep("'u0' must be a data.frame",
                      res[[1]]$message)) > 0)

u0 <- structure(list(S  = c(0, 1, 2, 3, 4, 5),
                     E  = c(0, 0, 0, 0, 0, 0),
                     I  = c(0, 0, 0, 0, 0, 0),
                     R  = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S", "E", "I", "R"),
                row.names = c(NA, -6L), class = "data.frame")

## Check missing columns in u0
res <- tools::assertError(SEIR(u0 = u0[, c("E", "I", "R"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SEIR(u0 = u0[, c("S", "I", "R"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SEIR(u0 = u0[, c("S", "E", "R"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(SEIR(u0 = u0[, c("S", "E", "I"), drop = FALSE]))
stopifnot(length(grep("Missing columns in u0",
                      res[[1]]$message)) > 0)

## Check missing beta
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               epsilon = 0.5,
                               gamma   = 0.5))
stopifnot(length(grep("'beta' is missing",
                      res[[1]]$message)) > 0)

## Check missing epsilon
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               gamma   = 0.5))
stopifnot(length(grep("'epsilon' is missing",
                      res[[1]]$message)) > 0)

## Check missing gamma
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = 0.5))
stopifnot(length(grep("'gamma' is missing",
                      res[[1]]$message)) > 0)

## Check non-numeric beta
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = "0.5",
                               epsilon = 0.3,
                               gamma   = 0.1))
stopifnot(length(grep("'beta' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric epsilon
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = "0.3",
                               gamma   = 0.1))
stopifnot(length(grep("'epsilon' must be numeric",
                      res[[1]]$message)) > 0)

## Check non-numeric gamma
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = 0.3,
                               gamma   = "0.1"))
stopifnot(length(grep("'gamma' must be numeric",
                      res[[1]]$message)) > 0)

## Check that length of beta equals 1
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = c(0.5, 0.5),
                               epsilon = 0.3,
                               gamma   = 0.1))
stopifnot(length(grep("'beta' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of epsilon equals 1
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = c(0.3, 0.3),
                               gamma   = 0.1))
stopifnot(length(grep("'epsilon' must be of length 1",
                      res[[1]]$message)) > 0)

## Check that length of gamma equals 1
res <- tools::assertError(SEIR(u0      = u0,
                               tspan   = seq_len(10) - 1,
                               events  = NULL,
                               beta    = 0.5,
                               epsilon = 0.3,
                               gamma   = c(0.1, 0.1)))
stopifnot(length(grep("'gamma' must be of length 1",
                      res[[1]]$message)) > 0)

## Check extraction of data from 'suscpetible', 'infected' and
## 'recovered' compartments
model <- SEIR(u0      = u0,
              tspan   = seq_len(10) - 1,
              events  = NULL,
              beta    = 0,
              epsilon = 0,
              gamma   = 0)

result <- run(model, threads = 1)

S_expected <- structure(c(0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L,
                          1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L,
                          2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L,
                          3L, 4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L, 0L, 1L, 2L, 3L,
                          4L, 5L, 0L, 1L, 2L, 3L, 4L, 5L),
                        .Dim = c(6L, 10L),
                        .Dimnames = list(c("S", "S", "S", "S", "S", "S"),
                                         c("0", "1", "2", "3", "4", "5",
                                           "6", "7", "8", "9")))

S_observed <- U(result, compartments = "S", as.is = TRUE)
stopifnot(identical(S_observed, S_expected))

I_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L),
                        .Dimnames = list(c("I", "I", "I", "I", "I", "I"),
                                         c("0", "1", "2", "3", "4", "5",
                                           "6", "7", "8", "9")))

I_observed <- U(result, compartments = "I", as.is = TRUE)
stopifnot(identical(I_observed, I_expected))

R_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(6L, 10L),
                        .Dimnames = list(c("R", "R", "R", "R", "R", "R"),
                                         c("0", "1", "2", "3", "4", "5",
                                           "6", "7", "8", "9")))

R_observed <- U(result, compartments = "R", as.is = TRUE)
stopifnot(identical(R_observed, R_expected))

R_expected <- structure(list(
    Node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L,
             3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L,
             1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L),
    Time = c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L,
             1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L,
             4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 7L,
             7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L,
             9L),
    R = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)),
    .Names = c("Node",  "Time", "R"),
    class = "data.frame", row.names = c(NA, -60L))
R_observed <- U(result, compartments = "R")
stopifnot(identical(R_observed, R_expected))

## Extract the number of recovered individuals in the first node after
## each time step in the simulation
R_expected <- structure(list(
    Node = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    Time = 0:9L,
    R = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)),
    .Names = c("Node", "Time", "R"),
    row.names = c(NA, -10L),
    class = "data.frame")
R_observed <- U(result, compartments = "R", i = 1)
stopifnot(identical(R_observed, R_expected))

R_expected <-structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                       .Dim = c(1L, 10L),
                       .Dimnames = list(c("R"),
                                        c("0", "1", "2", "3", "4", "5",
                                          "6", "7", "8", "9")))
R_observed <- U(result, compartments = "R", i = 1, as.is = TRUE)
stopifnot(identical(R_observed, R_expected))

## Extract the number of recovered individuals in the first and third
## node after each time step in the simulation
R_expected <- structure(list(
    Node = c(1L, 3L, 1L, 3L, 1L, 3L, 1L, 3L, 1L, 3L,
             1L, 3L, 1L, 3L, 1L, 3L, 1L, 3L, 1L, 3L),
    Time = c(0L, 0L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L,
             5L, 5L, 6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L),
    R = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)),
    .Names = c("Node", "Time", "R"),
    row.names = c(NA, -20L),
    class = "data.frame")
R_observed <- U(result, compartments = "R", i = c(1, 3))
stopifnot(identical(R_observed, R_expected))

R_expected <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                          0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                        .Dim = c(2L, 10L),
                        .Dimnames = list(c("R", "R"),
                                         c("0", "1", "2", "3", "4", "5",
                                           "6", "7", "8", "9")))
R_observed <- U(result, compartments = "R", i = c(1, 3), as.is = TRUE)
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
result <- run(model, threads = 1)

U_expected <- structure(list(
    Node = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L,
             3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L,
             1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
             5L, 6L),
    Time = c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L,
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
          649L)),
    .Names = c("Node", "Time", "S", "E", "I", "R"),
    row.names = c(NA, -60L),
    class = "data.frame")

U_observed <- U(result)
stopifnot(identical(U_observed, U_expected))
U_observed <- U(result, compartments = c("S", "E", "I", "R"), i = 1:6)
stopifnot(identical(U_observed, U_expected))

U_expected <- structure(list(
    Node = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    Time = 0:9, S = 110:119, I = 130:139),
    .Names = c("Node", "Time", "S", "I"),
    row.names = c(NA, -10L), class = "data.frame")
U_observed <- U(result, compartments = c("S", "I"), i = 1)
stopifnot(identical(U_observed, U_expected))

U_expected <-structure(list(
    Node = c(2L, 4L, 2L, 4L, 2L, 4L, 2L, 4L, 2L, 4L,
             2L, 4L, 2L, 4L, 2L, 4L, 2L, 4L, 2L, 4L),
    Time = c(0L, 0L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L,
             5L, 5L, 6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L),
    E = c(220L, 420L, 221L, 421L, 222L, 422L, 223L, 423L, 224L, 424L,
          225L, 425L, 226L, 426L, 227L, 427L, 228L, 428L, 229L, 429L),
    R = c(240L, 440L, 241L, 441L, 242L, 442L, 243L, 443L, 244L, 444L,
          245L, 445L, 246L, 446L, 247L, 447L, 248L, 448L, 249L, 449L)),
    .Names = c("Node", "Time", "E", "R"), row.names = c(NA, -20L),
    class = "data.frame")
U_observed <- U(result, compartments = c("E", "R"), i = c(2, 4))
stopifnot(identical(U_observed, U_expected))

U_expected <- structure(
    c(220L, 240L, 420L, 440L, 221L, 241L, 421L, 441L, 222L,
      242L, 422L, 442L, 223L, 243L, 423L, 443L, 224L, 244L, 424L, 444L,
      225L, 245L, 425L, 445L, 226L, 246L, 426L, 446L, 227L, 247L, 427L,
      447L, 228L, 248L, 428L, 448L, 229L, 249L, 429L, 449L),
    .Dim = c(4L, 10L),
    .Dimnames = list(c("E", "R", "E", "R"),
                     c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")))
U_observed <- U(result, compartments = c("E", "R"), i = c(2, 4), as.is = TRUE)
stopifnot(identical(U_observed, U_expected))

## Check 'V'
res <- tools::assertError(V(result))
stopifnot(length(grep("No continuous variables defined in 'model'",
                      res[[1]]$message)) > 0)

## Check data
stopifnot(identical(nrow(events_SEIR()), 466692L))
stopifnot(identical(nrow(u0_SEIR()), 1600L))

## Check SEIR plot method
pdf_file <- tempfile(fileext = ".pdf")
pdf(pdf_file)
plot(result)
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
res <- tools::assertError(.Call("SEIR_run", NULL, NULL, NULL, NULL,
                                PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid model.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(.Call("SEIR_run", "SEIR", NULL, NULL, NULL,
                                PACKAGE = "SimInf"))
stopifnot(length(grep("Invalid model.",
                      res[[1]]$message)) > 0)
