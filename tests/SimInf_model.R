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

## Specify the number of threads to use.
set_num_threads(1)

## For debugging
sessionInfo()

## Initialize test data
S <- SimInf:::init_sparse_matrix(Matrix(c(-1,  0,  0,
                                          1,  0,  0,
                                          0, -1,  0,
                                          0,  1,  0,
                                          0,  0, -1,
                                          0,  0,  1),
                                        nrow   = 6,
                                        ncol   = 3,
                                        byrow  = TRUE,
                                        dimnames = list(LETTERS[1:6], NULL)))

Nn <- 6L

G <- SimInf:::init_sparse_matrix(Matrix(c(1, 0, 0,
                                          0, 1, 0,
                                          0, 0, 1),
                                        nrow   = 3,
                                        ncol   = 3,
                                        byrow  = TRUE,
                                        dimnames = list(
                                            c("A -> B", "C -> D", "E -> F"),
                                            NULL)))

u0 <- structure(c(0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2,
                  0, 3, 0, 3, 0, 3, 0, 4, 0, 4, 0, 4, 0, 5, 0, 5, 0, 5, 0),
                .Dim = c(6L, 6L))
storage.mode(u0) <- "integer"

U <- matrix(nrow = 0, ncol = 0)
storage.mode(U) <- "integer"

## Check valid_SimInf_model_object
m <- SISe(u0 = c(S = 10, I = 0), tspan = 1:10, phi = 0,
          upsilon = 0.1, gamma = 0.1, alpha = 1.0, beta_t1 = 0.1,
          beta_t2 = 0.1, beta_t3 = 0.1, beta_t4 = 0.1, end_t1  = 91,
          end_t2  = 182, end_t3  = 273, end_t4  = 365, epsilon = 0.1)
stopifnot(isTRUE(SimInf:::valid_SimInf_model_object(m)))

## Check valid_SimInf_model_object with invalid tspan.
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
m@tspan <- integer(0)
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "Input time-span must be a double vector."))
m@tspan <- numeric(0)
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "Input time-span must be an increasing vector."))

## Check valid_SimInf_model_object with invalid u0.
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
storage.mode(m@u0) <- "double"
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "Initial state 'u0' must be an integer matrix."))
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
m@u0[1, 1] <- -1L
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "Initial state 'u0' has negative elements."))

## Check valid_SimInf_model_object with invalid U.
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
storage.mode(m@U) <- "double"
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "Output state 'U' must be an integer matrix."))
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
m@U <- matrix(-1L)
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "Output state 'U' has negative elements."))

## Check valid_SimInf_model_object with invalid v0.
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
m@v0 <- matrix(1L)
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "Initial model state 'v0' must be a double matrix."))
m@v0 <- matrix(1)
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "'v0' must have rownames."))
m@v0 <- matrix(c(1, 1), ncol = 2)
rownames(m@v0) <- "test"
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "The number of nodes in 'u0' and 'v0' must match."))

## Check valid_SimInf_model_object with invalid V.
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
storage.mode(m@V) <- "integer"
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "Output model state 'V' must be a double matrix."))

## Check valid_SimInf_model_object with invalid S.
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
m@S@x <- m@S@x * 0.5
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "'S' matrix must be an integer matrix."))

## Check valid_SimInf_model_object with different rownames in S and E.
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1,
         events = data.frame(event = 1, node = 1, n = 1, time = 1, dest = 0,
                             proportion = 0, select = 1, shift = 0))
rownames(m@events@E) <- NULL
stopifnot(identical(
    SimInf:::valid_SimInf_model_object(m),
    "'S' and 'E' must have rownames matching the compartments."))
rownames(m@events@E) <- rownames(m@S)
rownames(m@events@E)[1] <- "Z"
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "'S' and 'E' must have identical compartments."))

## Check valid_SimInf_model_object with invalid G.
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
m@G <- m@G[1, 1:2, drop = FALSE]
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "Wrong size of dependency graph."))

m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
rownames(m@G) <- NULL
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "'G' must have rownames that specify transitions."))

m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
rownames(m@G)[1] <- ""
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "'G' must have rownames that specify transitions."))

m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
rownames(m@G)[1] <- "A"
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "'G' rownames have invalid transitions."))

m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
rownames(m@G)[1] <- "Z -> beta*Z*I/(S+I+R) -> I"
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "'G' and 'S' must have identical compartments."))

## Check valid_SimInf_model_object with invalid ldata.
m <- SIR(u0 = data.frame(S = 10, I = 0, R = 0),
         tspan = 1:10, beta = 0.1, gamma = 0.1)
m@ldata <- matrix(1L)
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "'ldata' matrix must be a double matrix."))
m@ldata <- matrix(c(1, 2), ncol = 2)
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "The number of nodes in 'u0' and 'ldata' must match."))

## Check valid_SimInf_model_object with invalid gdata.
m <- SISe(u0 = data.frame(S = 10, I = 0), tspan = 1:10, phi = 0,
          upsilon = 0.1, gamma = 0.1, alpha = 1.0, beta_t1 = 0.1,
          beta_t2 = 0.1, beta_t3 = 0.1, beta_t4 = 0.1, end_t1  = 91,
          end_t2  = 182, end_t3  = 273, end_t4  = 365, epsilon = 0.1)
storage.mode(m@gdata) <- "integer"
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "'gdata' must be a double vector."))

## Check v0
m <- SISe(u0 = data.frame(S = 10, I = 0), tspan = 1:10, phi = 0,
          upsilon = 0.1, gamma = 0.1, alpha = 1.0, beta_t1 = 0.1,
          beta_t2 = 0.1, beta_t3 = 0.1, beta_t4 = 0.1, end_t1  = 91,
          end_t2  = 182, end_t3  = 273, end_t4  = 365, epsilon = 0.1)

storage.mode(m@v0) <- "integer"
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "Initial model state 'v0' must be a double matrix."))

storage.mode(m@v0) <- "double"
rownames(m@v0) <- NULL
stopifnot(identical(SimInf:::valid_SimInf_model_object(m),
                    "'v0' must have rownames."))

res <- SimInf_model(G     = G,
                    S     = S,
                    U     = U,
                    ldata = matrix(rep(0, Nn), nrow = 1),
                    tspan = c(1, 2),
                    u0    = u0,
                    v0    = matrix(rep(1L, Nn), nrow = 1,
                                   dimnames = list("phi")))

## Check tspan
res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = numeric(0),
                       u0    = u0))
check_error(res, "Input time-span must be an increasing vector.", FALSE)

res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = as.numeric(c(3, 2, 1)),
                       u0    = u0))
check_error(res, "Input time-span must be an increasing vector.", FALSE)

res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = c(1, NA, 3),
                       u0    = u0))
check_error(res, "Input time-span must be an increasing vector.", FALSE)

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
res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = as.numeric(1:10),
                       u0    = u0 * -1L))
check_error(res, "Initial state 'u0' has negative elements.", FALSE)

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
res <- assertError(SimInf_model(G     = G,
                                S     = S,
                                U     = U,
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                tspan = as.numeric(1:10),
                                u0    = u0_double))
check_error(res, "'u0' must be an integer matrix.")

## Change u0 to vector. Should raise an error.
res <- assertError(SimInf_model(G     = G,
                                S     = S,
                                U     = U,
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                tspan = as.numeric(1:10),
                                u0    = as.numeric(u0)))
check_error(res, "The number of rows in 'u0' and 'S' must match.", FALSE)
check_error(res, "The number of nodes in 'u0' and 'ldata' must match.", FALSE)

## Check S
res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S * 1.1,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = as.numeric(1:10),
                       u0    = u0))
check_error(res, "'S' matrix must be an integer matrix.", FALSE)

## Check G
## Error: Wrong size of dependency graph
res <- assertError(new("SimInf_model",
                       G     = G[-1, ],
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = as.numeric(1:10),
                       u0    = u0))
check_error(res, "Wrong size of dependency graph.", FALSE)

## Check specication of transition
rownames(G) <- NULL
res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = as.numeric(1:10),
                       u0    = u0))
check_error(res, "'G' must have rownames that specify transitions.", FALSE)

rownames(G) <- c("", "  ", "E -> F")
res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = as.numeric(1:10),
                       u0    = u0))
check_error(res, "'G' must have rownames that specify transitions.", FALSE)

rownames(G) <- c("A -> B", "C -> D", "E ->")
res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = as.numeric(1:10),
                       u0    = u0))
check_error(res, "'G' rownames have invalid transitions.", FALSE)

rownames(G) <- c("A -> B", "C -> D", "E ->")
res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = as.numeric(1:10),
                       u0    = u0))
check_error(res, "'G' rownames have invalid transitions.", FALSE)

rownames(G) <- c("A -> B", "C -> D", "E -> G")
res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       tspan = as.numeric(1:10),
                       u0    = u0))
check_error(res, "'G' and 'S' must have identical compartments", FALSE)
rownames(G) <- c("A -> B", "C -> D", "E -> F")

## Check gdata
res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = matrix(rep(0, Nn), nrow = 1),
                       gdata = 1L,
                       tspan = as.numeric(1:10),
                       u0    = u0))
check_error(res, "'gdata' must be a double vector.", FALSE)

## Check ldata
ldata <- matrix(rep(0, Nn), nrow = 1)
ldata <- ldata[, 1:3, drop = FALSE]

## Wrong size of ldata matrix
res <- assertError(new("SimInf_model",
                       G     = G,
                       S     = S,
                       U     = U,
                       ldata = ldata,
                       tspan = as.numeric(1:10),
                       u0    = u0))
check_error(res, "The number of nodes in 'u0' and 'ldata' must match.", FALSE)

## Check initial state
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

## 'u0' is NULL
res <- assertError(SimInf_model())
check_error(res, "'u0' must be an integer matrix.")

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

stopifnot(identical(show(model), model))
show_observed <- capture.output(show(model))

stopifnot(identical(show_observed[1:5], show_expected))

## Check summary method without events
summary(run(model))

## Check first lines of show method with events
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

events <- data.frame(
    event      = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
    time       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    node       = c(2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    dest       = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    n          = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    proportion = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    select     = c(4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6),
    shift      = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

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

show_observed <- capture.output(show(model))

stopifnot(identical(
    show_observed[1:5],
    c("Model: SISe3",
      "Number of nodes: 6",
      "Number of transitions: 6",
      "Number of scheduled events: 15",
      "")))

stopifnot(identical(
    show_observed[22:28],
    c("Local data",
      "----------",
      " Parameter Value",
      " end_t1     91  ",
      " end_t2    182  ",
      " end_t3    273  ",
      " end_t4    365  ")))

## Check summary method with events
summary(run(model))

## Check U. Change storage mode of U to double.
## Should not raise error
U <- structure(c(
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L, 2L, 0L, 1L, 1L,
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

res <- assertError(SimInf_model(G     = G,
                                S     = S,
                                U     = U_double,
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                tspan = as.numeric(1:10),
                                u0    = u0))
check_error(res, "'U' must be an integer matrix.")

## Check U. Should not raise an error if U is an integer vector of length 0
SimInf_model(G     = G,
             S     = S,
             U     = integer(0),
             ldata = matrix(rep(0, Nn), nrow = 1),
             tspan = as.numeric(1:10),
             u0    = u0)

## Check U. Should raise error if U is an integer vector of length > 0
res <- assertError(SimInf_model(G     = G,
                                S     = S,
                                U     = c(1L),
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                tspan = as.numeric(1:10),
                                u0    = u0))
check_error(res, "'U' must be equal to a 0 x 0 matrix.")

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

res <- assertError(SimInf_model(G     = G,
                                S     = S,
                                U     = U,
                                V     = V_character,
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                tspan = as.numeric(1:10),
                                u0    = u0))
check_error(res, "'V' must be a double matrix.")

## Check V. Should raise error if V is a vector of length > 0
res <- assertError(SimInf_model(G     = G,
                                S     = S,
                                U     = U,
                                V     = 1,
                                ldata = matrix(rep(0, Nn), nrow = 1),
                                tspan = as.numeric(1:10),
                                u0    = u0))
check_error(res, "'V' must be equal to a 0 x 0 matrix.")

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
res <- assertError(plot(model))
check_error(res, "Please run the model first, the trajectory is empty.")

## Check that the SimInf_model initialisation fails if the events
## argument is not either NULL or a data.frame
u0 <- data.frame(S_1 = c(0, 1, 2, 3, 4, 5),
                 I_1 = c(0, 0, 0, 0, 0, 0),
                 S_2 = c(0, 1, 2, 3, 4, 5),
                 I_2 = c(0, 0, 0, 0, 0, 0),
                 S_3 = c(0, 1, 2, 3, 4, 5),
                 I_3 = c(0, 0, 0, 0, 0, 0))

res <- assertError(model <- SISe3(u0        = u0,
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
check_error(res, "'events' must be NULL or a data.frame.")

## Check arguments to 'trajectory' method
u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
result <- run(model)
res <- assertError(trajectory(result, compartments = c("A", "S")))
check_error(res, "Non-existing compartment(s) in model: 'A'.")
res <- assertError(trajectory(result, index = c("A", "S")))
check_error(
    res,
    "The node index must be an integer > 0 and <= number of nodes.")
res <- assertError(trajectory(result, index = 3.4))
check_error(
    res,
    "The node index must be an integer > 0 and <= number of nodes.")
res <- assertError(trajectory(result, index = 0))
check_error(
    res,
    "The node index must be an integer > 0 and <= number of nodes.")
res <- assertError(trajectory(result, index = 10))
check_error(
    res,
    "The node index must be an integer > 0 and <= number of nodes.")

## Check arguments to 'prevalence' method
u0 <- data.frame(S = 100:105, I = 1:6, R = rep(0, 6))
model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
result <- run(model)
res <- assertError(prevalence(result, A + S ~ S))
check_error(res, "Non-existing compartment(s) in model: 'A'.")
res <- assertError(prevalence(result, S ~ A + S))
check_error(res, "Non-existing compartment(s) in model: 'A'.")
res <- assertError(prevalence(result, I ~ S + I + R, i = c("A", "S")))
check_error(
    res,
    "The node index must be an integer > 0 and <= number of nodes.")
res <- assertError(prevalence(result, I ~ S + I + R, i = 3.4))
check_error(
    res,
    "The node index must be an integer > 0 and <= number of nodes.")
res <- assertError(prevalence(result, I ~ S + I + R, i = 0))
check_error(
    res,
    "The node index must be an integer > 0 and <= number of nodes.")
res <- assertError(prevalence(result, I ~ S + I + R, i = 10))
check_error(
    res,
    "The node index must be an integer > 0 and <= number of nodes.")

## Check 'gdata'
model <- SISe(u0 = data.frame(S = 10, I = 0), tspan = 1:10, phi = 0,
              upsilon = 1, gamma = 2, alpha = 3, beta_t1 = 4,
              beta_t2 = 5, beta_t3 = 6, beta_t4 = 7, end_t1  = 91,
              end_t2  = 182, end_t3  = 273, end_t4  = 365, epsilon = 8)

stopifnot(identical(gdata(model),
                    c(upsilon = 1, gamma = 2, alpha = 3, beta_t1 = 4,
                      beta_t2 = 5, beta_t3 = 6, beta_t4 = 7, epsilon = 8)))

gdata(model, "epsilon") <- 9

stopifnot(identical(gdata(model),
                    c(upsilon = 1, gamma = 2, alpha = 3, beta_t1 = 4,
                      beta_t2 = 5, beta_t3 = 6, beta_t4 = 7, epsilon = 9)))

res <- assertError(gdata(model) <- 6)
check_error(res, "Missing 'parameter' argument.")

res <- assertError(gdata(model, "epsilon") <- "6")
check_error(res, "'value' argument must be a numeric.")

res <- assertError(gdata(model, 5) <- 6)
check_error(res, "'parameter' argument must be a character.")

res <- assertError("gdata<-" (model, "epsilon"))
check_error(res, "Missing 'value' argument.")

## Check 'ldata'
model@ldata <- matrix(1, dimnames = list("test", NULL))
res <- assertError(ldata(model))
check_error(res, "Missing 'node' argument.")

res <- assertError(ldata(model, "0"))
check_error(res, "Invalid 'node' argument.")

res <- assertError(ldata(model, 0))
check_error(res, "Invalid 'node' argument.")

res <- assertError(ldata(model, c(0, 0)))
check_error(res, "Invalid 'node' argument.")

ldata(model, 1)
