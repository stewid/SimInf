## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
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

library("SimInf")
library("Matrix")

## For debugging
sessionInfo()

## Initialize test data
S <- Matrix(c(-1,  0,  0,
               1,  0,  0,
               0, -1,  0,
               0,  1,  0,
               0,  0, -1,
               0,  0,  1),
            nrow   = 6,
            ncol   = 3,
            byrow  = TRUE,
            sparse = TRUE)
rownames(S) <- LETTERS[1:6]

Nn <- 6L

G <- as(Matrix(c(1, 0, 0,
                 0, 1, 0,
                 0, 0, 1),
               nrow   = 3,
               ncol   = 3,
               byrow  = TRUE,
               sparse = TRUE),
        "dgCMatrix")
rownames(G) <- c("A -> B", "C -> D", "E -> F")

u0 <- structure(c(0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2,
                  0, 3, 0, 3, 0, 3, 0, 4, 0, 4, 0, 4, 0, 5, 0, 5, 0, 5, 0),
                .Dim = c(6L, 6L))
storage.mode(u0) <- "integer"

U <- matrix(nrow = 0, ncol = 0)
storage.mode(U) <- "integer"

## Check tspan
res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = as.numeric(1),
                              u0    = u0))
stopifnot(length(grep("Input time-span must be an increasing vector.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = as.numeric(c(3, 2, 1)),
                              u0    = u0))
stopifnot(length(grep("Input time-span must be an increasing vector.",
                      res[[1]]$message)) > 0)

res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = c(1, NA, 3),
                              u0    = u0))
stopifnot(length(grep("Input time-span must be an increasing vector.",
                      res[[1]]$message)) > 0)

## Check that tspan can be a Date vector
res <- SimInf_model(G     = G,
                    S     = S,
                    U     = U,
                    ldata = matrix(rep(0, Nn), nrow = 1),
                    tspan = as.Date(c("2017-01-01",
                                      "2017-01-02",
                                      "2017-01-03")),
                    u0    = u0)
stopifnot(identical(res@tspan,
                    structure(c(1, 2, 3),
                              .Names = c("2017-01-01",
                                         "2017-01-02",
                                         "2017-01-03"))))

## Check u0
res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = as.numeric(1:10),
                              u0    = u0 * -1L))
stopifnot(length(grep("Initial state 'u0' has negative elements.",
                      res[[1]]$message)) > 0)

## Change storage mode of u0 to double.
## Should not raise error
u0_double <- u0
storage.mode(u0_double) <- "double"
SimInf_model(G     = G,
             S     = S,
             U     = U,
             ldata = matrix(rep(0, Nn), nrow = 1),
             tspan = as.numeric(1:10),
             u0    = u0_double)

## Change storage mode of u0 to double and change to non-integer values.
## Should raise error
u0_double <- u0
storage.mode(u0_double) <- "double"
u0_double <- 1.2 * u0_double
res <- tools::assertError(SimInf_model(G     = G,
                                       S     = S,
                                       U     = U,
                                       ldata = matrix(rep(0, Nn), nrow = 1),
                                       tspan = as.numeric(1:10),
                                       u0    = u0_double))
stopifnot(length(grep("u0 must be an integer matrix",
                      res[[1]]$message)) > 0)

## Check S
res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S * 1.1,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = as.numeric(1:10),
                              u0    = u0))
stopifnot(length(grep("'S' matrix must be an integer matrix.",
                      res[[1]]$message)) > 0)

## Check G
## Error: Wrong size of dependency graph
res <- tools::assertError(new("SimInf_model",
                              G     = G[-1,],
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = as.numeric(1:10),
                              u0    = u0))
stopifnot(length(grep("Wrong size of dependency graph.",
                      res[[1]]$message)) > 0)

## Check specication of transition
rownames(G) <- NULL
res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = as.numeric(1:10),
                              u0    = u0))
stopifnot(length(grep("'G' must have rownames that specify transitions.",
                      res[[1]]$message)) > 0)

rownames(G) <- c("", "  ", "E -> F")
res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = as.numeric(1:10),
                              u0    = u0))
stopifnot(length(grep("'G' must have rownames that specify transitions.",
                      res[[1]]$message)) > 0)

rownames(G) <- c("A -> B", "C -> D", "E ->")
res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = as.numeric(1:10),
                              u0    = u0))
stopifnot(length(grep("'G' rownames have invalid transitions.",
                      res[[1]]$message)) > 0)

rownames(G) <- c("A -> B", "C -> D", "E ->")
res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = as.numeric(1:10),
                              u0    = u0))
stopifnot(length(grep("'G' rownames have invalid transitions.",
                      res[[1]]$message)) > 0)

rownames(G) <- c("A -> B", "C -> D", "E -> G")
res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              tspan = as.numeric(1:10),
                              u0    = u0))
stopifnot(length(grep("'G' and 'S' must have identical compartments",
                      res[[1]]$message)) > 0)
rownames(G) <- c("A -> B", "C -> D", "E -> F")

## Check gdata
res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = matrix(rep(0, Nn), nrow = 1),
                              gdata = 1L,
                              tspan = as.numeric(1:10),
                              u0    = u0))
stopifnot(length(grep("'gdata' must be a double vector.",
                      res[[1]]$message)) > 0)

## Check ldata
ldata <- matrix(rep(0, Nn), nrow = 1)
ldata <- ldata[, 1:3, drop = FALSE]

## Wrong size of ldata matrix
res <- tools::assertError(new("SimInf_model",
                              G     = G,
                              S     = S,
                              U     = U,
                              ldata = ldata,
                              tspan = as.numeric(1:10),
                              u0    = u0))
stopifnot(length(grep("Wrong size of 'ldata' matrix.",
                      res[[1]]$message)) > 0)

## Check initial state
u0 <- structure(list(S_1 = c(0, 1, 2, 3, 4, 5),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 1, 2, 3, 4, 5),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 1, 2, 3, 4, 5),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L), class = "data.frame")

## 'u0' is NULL
res <- tools::assertError(SimInf_model())
stopifnot(length(grep("'u0' is NULL",
                      res[[1]]$message)) > 0)

## Check first lines of show method without events
model <- SISe(u0      = data.frame(S = 99, I = 1),
              tspan   = seq_len(1000) - 1,
              events  = NULL,
              phi     = 1,
              upsilon = 0.017,
              gamma   = 0.1,
              alpha   = 1,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0.000011)

show_expected <- c("Model: SISe",
                   "Number of nodes: 1",
                   "Number of transitions: 2",
                   "Number of scheduled events: 0",
                   "")

show_observed <- capture.output(show(model))

stopifnot(identical(show_observed[1:5], show_expected))

## Check summary method without events
summary(run(model))

## Check first lines of show method with events
u0 <- structure(list(S_1 = c(0, 1, 2, 3, 4, 5),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 1, 2, 3, 4, 5),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 1, 2, 3, 4, 5),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L),
                class = "data.frame")

events <- structure(list(
    event      = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node       = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6),
    shift      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
                    .Names = c("event", "time", "node", "dest",
                        "n", "proportion", "select", "shift"),
                    row.names = c(NA, -15L), class = "data.frame")

model <- SISe3(u0        = u0,
               tspan     = 0:10,
               events    = events,
               phi       = rep(1, 6),
               upsilon_1 = 1,
               upsilon_2 = 1,
               upsilon_3 = 1,
               gamma_1   = 1,
               gamma_2   = 1,
               gamma_3   = 1,
               alpha     = 1,
               beta_t1   = 1,
               beta_t2   = 1,
               beta_t3   = 1,
               beta_t4   = 1,
               end_t1    = 91,
               end_t2    = 182,
               end_t3    = 273,
               end_t4    = 365,
               epsilon   = 1)

show_expected <- c("Model: SISe3",
                   "Number of nodes: 6",
                   "Number of transitions: 6",
                   "Number of scheduled events: 15",
                   "")

show_observed <- capture.output(show(model))

stopifnot(identical(show_observed[1:5], show_expected))

## Check summary method with events
summary(run(model))

## Check U. Change storage mode of U to double.
## Should not raise error
U <- structure(c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L, 2L, 0L, 1L, 1L,
                 1L, 1L, 2L, 1L, 3L, 0L, 2L, 1L, 2L, 2L, 0L, 4L, 1L, 3L, 2L, 3L,
                 3L, 2L, 1L, 4L, 6L, 9L, 7L, 8L, 4L, 11L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 7L, 8L, 7L, 8L, 5L, 10L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 8L, 7L, 6L, 9L,
                 8L, 7L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 5L, 10L, 4L, 11L, 6L, 9L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 7L, 8L, 5L, 10L, 7L, 8L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L, 11L, 5L, 10L, 3L, 12L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 8L, 7L, 5L,
                 10L, 4L, 11L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 6L, 9L, 2L, 13L, 4L, 11L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 5L, 10L, 2L, 13L, 7L, 8L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 9L, 6L, 2L, 13L, 6L, 9L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
               .Dim = c(36L, 11L))

U_double <- U
storage.mode(U_double) <- "double"

SimInf_model(G     = G,
             S     = S,
             U     = U_double,
             ldata = matrix(rep(0, Nn), nrow = 1),
             tspan = as.numeric(1:10),
             u0    = u0)

## Check U. Change storage mode of U to double and change to non-integer values.
## Should raise error
U_double <- U
storage.mode(U_double) <- "double"
U_double <- U_double * 1.2

res <- tools::assertError(SimInf_model(G     = G,
                                       S     = S,
                                       U     = U_double,
                                       ldata = matrix(rep(0, Nn), nrow = 1),
                                       tspan = as.numeric(1:10),
                                       u0    = u0))
stopifnot(length(grep("U must be an integer",
                      res[[1]]$message)) > 0)

## Check U. Should not raise an error if U is an integer vector of length 0
SimInf_model(G     = G,
             S     = S,
             U     = integer(0),
             ldata = matrix(rep(0, Nn), nrow = 1),
             tspan = as.numeric(1:10),
             u0    = u0)

## Check U. Should raise error if U is an integer vector of length > 0
res <- tools::assertError(SimInf_model(G     = G,
                                       S     = S,
                                       U     = c(1L),
                                       ldata = matrix(rep(0, Nn), nrow = 1),
                                       tspan = as.numeric(1:10),
                                       u0    = u0))
stopifnot(length(grep("U must be equal to 0 x 0 matrix",
                      res[[1]]$message)) > 0)

## Check V. Change storage mode of V to double.
## Should not raise error
V <- structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1),
               .Dim = c(6L, 10L))

V_integer <- V
storage.mode(V) <- "integer"

SimInf_model(G     = G,
             S     = S,
             U     = U,
             V     = V_integer,
             ldata = matrix(rep(0, Nn), nrow = 1),
             tspan = as.numeric(1:10),
             u0    = u0)

## Check V. Change storage mode of V to character
## Should raise error
V_character <- V
storage.mode(V_character) <- "character"

res <- tools::assertError(SimInf_model(G     = G,
                                       S     = S,
                                       U     = U,
                                       V     = V_character,
                                       ldata = matrix(rep(0, Nn), nrow = 1),
                                       tspan = as.numeric(1:10),
                                       u0    = u0))
stopifnot(length(grep("V must be numeric",
                      res[[1]]$message)) > 0)

## Check V. Should raise error if V is a vector of length > 0
res <- tools::assertError(SimInf_model(G     = G,
                                       S     = S,
                                       U     = U,
                                       V     = 1,
                                       ldata = matrix(rep(0, Nn), nrow = 1),
                                       tspan = as.numeric(1:10),
                                       u0    = u0))
stopifnot(length(grep("V must be equal to 0 x 0 matrix",
                      res[[1]]$message)) > 0)

## Check V. Should not raise an error if V is an integer vector of length 0
SimInf_model(G     = G,
             S     = S,
             U     = U,
             V     = integer(0),
             ldata = matrix(rep(0, Nn), nrow = 1),
             tspan = as.numeric(1:10),
             u0    = u0)

## Check that plot raises an error if the model hasn't run
model <- SISe(u0      = data.frame(S = 99, I = 1),
              tspan   = seq_len(1000) - 1,
              events  = NULL,
              phi     = 1,
              upsilon = 0.017,
              gamma   = 0.1,
              alpha   = 1,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0.000011)
res <- tools::assertError(plot(model))
stopifnot(length(grep("Please run the model first, the 'U' matrix is empty",
                      res[[1]]$message)) > 0)

## Check that the SimInf_model initialisation fails if the events
## argument is not either NULL or a data.frame
u0 <- structure(list(S_1 = c(0, 1, 2, 3, 4, 5),
                     I_1 = c(0, 0, 0, 0, 0, 0),
                     S_2 = c(0, 1, 2, 3, 4, 5),
                     I_2 = c(0, 0, 0, 0, 0, 0),
                     S_3 = c(0, 1, 2, 3, 4, 5),
                     I_3 = c(0, 0, 0, 0, 0, 0)),
                .Names = c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3"),
                row.names = c(NA, -6L),
                class = "data.frame")

res <- tools::assertError(model <- SISe3(u0        = u0,
                                         tspan     = 0:10,
                                         events    = "events",
                                         phi       = rep(1, 6),
                                         upsilon_1 = 1,
                                         upsilon_2 = 1,
                                         upsilon_3 = 1,
                                         gamma_1   = 1,
                                         gamma_2   = 1,
                                         gamma_3   = 1,
                                         alpha     = 1,
                                         beta_t1   = 1,
                                         beta_t2   = 1,
                                         beta_t3   = 1,
                                         beta_t4   = 1,
                                         end_t1    = 91,
                                         end_t2    = 182,
                                         end_t3    = 273,
                                         end_t4    = 365,
                                         epsilon   = 1))
stopifnot(length(grep("'events' must be NULL or a data.frame",
                      res[[1]]$message)) > 0)

## Check that tspan names are copied to U and V colnames
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:10,
             beta = 0.16,
             gamma = 0.077)
result <- run(model, threads = 1)
stopifnot(identical(colnames(trajectory(result, as.is = TRUE)), as.character(1:10)))
stopifnot(identical(colnames(result@V), as.character(1:10)))

tspan <- seq(as.Date("2016-01-01"), as.Date("2016-01-10"), by = 1)
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = tspan,
             beta = 0.16,
             gamma = 0.077)
result <- run(model, threads = 1)
stopifnot(identical(colnames(trajectory(result, as.is = TRUE)), as.character(tspan)))
stopifnot(identical(colnames(result@V), as.character(tspan)))

u0 <- data.frame(S = 100:105, I = 1:6)
model <- SISe(u0 = u0, tspan = 1:10,
              phi = rep(0, 6),
              upsilon = 0.017,
              gamma   = 0.1,
              alpha   = 1,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0.000011)

U(model) <- data.frame(time = c(5, 6, 7, 8, 9, 10),
                       node = c(1, 1, 2, 2, 3, 3),
                       S = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
                       I = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE))
V(model) <- data.frame(time = 5:10,
                       node = 1:6,
                       V1 = rep(TRUE, 6))
result <- run(model, threads = 1)
stopifnot(identical(colnames(trajectory(result, as.is = TRUE)), as.character(1:10)))
stopifnot(identical(colnames(trajectory(result, "V1", as.is = TRUE)), as.character(1:10)))

tspan <- seq(as.Date("2016-01-01"), as.Date("2016-01-10"), by = 1)
u0 <- data.frame(S = 100:105, I = 1:6)
model <- SISe(u0 = u0, tspan = tspan,
              phi = rep(0, 6),
              upsilon = 0.017,
              gamma   = 0.1,
              alpha   = 1,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0.000011)

U(model) <- data.frame(time = c(5, 6, 7, 8, 9, 10),
                       node = c(1, 1, 2, 2, 3, 3),
                       S = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
                       I = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE))
V(model) <- data.frame(time = 5:10,
                       node = 1:6,
                       V1 = rep(TRUE, 6))
V(model) <- data.frame(time = 5:10,
                       node = 1:6,
                       V1 = rep(TRUE, 6))
result <- run(model, threads = 1)
stopifnot(identical(colnames(trajectory(result, as.is = TRUE)), as.character(tspan)))
stopifnot(identical(colnames(trajectory(result, "V1", as.is = TRUE)), as.character(tspan)))

## Check arguments to 'U' method
u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
result <- run(model, threads = 1)
res <- tools::assertError(trajectory(result, compartments = c("A", "S")))
stopifnot(length(grep("Non-existing compartment[(]s[)] in model: 'A'",
                      res[[1]]$message)) > 0)
res <- tools::assertError(trajectory(result, i = c("A", "S")))
stopifnot(length(grep("'i' must be integer",
                      res[[1]]$message)) > 0)
res <- tools::assertError(trajectory(result, i = 3.4))
stopifnot(length(grep("'i' must be integer",
                      res[[1]]$message)) > 0)
res <- tools::assertError(trajectory(result, i = 0))
stopifnot(length(grep("'i' must be integer > 0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(trajectory(result, i = 10))
stopifnot(length(grep("'i' must be integer <= number of nodes",
                      res[[1]]$message)) > 0)

## Check arguments to 'prevalence' method
u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
result <- run(model, threads = 1)
res <- tools::assertError(prevalence(result, A+S~S))
stopifnot(length(grep("Non-existing compartment[(]s[)] in model: 'A'",
                      res[[1]]$message)) > 0)
res <- tools::assertError(prevalence(result, S~A+S))
stopifnot(length(grep("Non-existing compartment[(]s[)] in model: 'A'",
                      res[[1]]$message)) > 0)
res <- tools::assertError(prevalence(result, I~S+I+R, i = c("A", "S")))
stopifnot(length(grep("'i' must be integer",
                      res[[1]]$message)) > 0)
res <- tools::assertError(prevalence(result, I~S+I+R, i = 3.4))
stopifnot(length(grep("'i' must be integer",
                      res[[1]]$message)) > 0)
res <- tools::assertError(prevalence(result, I~S+I+R, i = 0))
stopifnot(length(grep("'i' must be integer > 0",
                      res[[1]]$message)) > 0)
res <- tools::assertError(prevalence(result, I~S+I+R, i = 10))
stopifnot(length(grep("'i' must be integer <= number of nodes",
                      res[[1]]$message)) > 0)

## Check 'gdata'
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 2, gamma = 4)

stopifnot(identical(gdata(model), structure(c(2, 4), .Names = c("beta", "gamma"))))
gdata(model, "beta") <- 6
stopifnot(identical(gdata(model), structure(c(6, 4), .Names = c("beta", "gamma"))))
tools::assertError(gdata(model) <- 6)
tools::assertError(gdata(model, "beta") <- "6")
