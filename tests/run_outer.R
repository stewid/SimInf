## SimInf, a framework for stochastic disease spread simulations
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
source("util/check.R")

## For debugging
sessionInfo()

## Define some variables
tol <- 1e-8
x <- seq(from = 0.95, to = 1.05, by = 0.01)
y <- seq(from = 0.95, to = 1.05, by = 0.01)

## Check gdata names
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
names(model@gdata) <- NULL
res <- tools::assertError(
    run_outer(x, y, model, alpha ~ upsilon, function(m) 1))
check_error(res, "'names(model@gdata)' is NULL.")

## Check formula argument
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
res <- tools::assertError(run_outer(x, y, model, NULL, function(m) 1))
check_error(res, "'formula' argument is NULL.")

## Check FUN argument
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
res <- tools::assertError(run_outer(x, y, model, a ~ b, NULL))
check_error(res, "'FUN' argument is NULL.")

## Check lhs
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
res <- tools::assertError(
    run_outer(x, y, model,  ~ upsilon, function(m) 1))
check_error(res, "Invalid parameters on the left side of the formula.")

res <- tools::assertError(
    run_outer(x, y, model, dummy ~ upsilon, function(m) 1))
check_error(res, "Unmatched parameters on the left hand side of the formula.")

## Check rhs
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
res <- tools::assertError(
    run_outer(x, y, model, alpha ~ upsilon:alpha, function(m) 1))
check_error(res, "Invalid parameters on the right side of the formula.")

res <- tools::assertError(
    run_outer(x, y, model, alpha ~ dummy, function(m) 1))
check_error(res, "Unmatched parameters on the right side of the formula.")

## Check run_outer
z_exp <- structure(
    c(0.008745225, 0.00883728, 0.008929335, 0.00902139,
      0.009113445, 0.0092055, 0.009297555, 0.00938961, 0.009481665,
      0.00957372, 0.009665775, 0.00883728, 0.008930304, 0.009023328,
      0.009116352, 0.009209376, 0.0093024, 0.009395424, 0.009488448,
      0.009581472, 0.009674496, 0.00976752, 0.008929335, 0.009023328,
      0.009117321, 0.009211314, 0.009305307, 0.0093993, 0.009493293,
      0.009587286, 0.009681279, 0.009775272, 0.009869265, 0.00902139,
      0.009116352, 0.009211314, 0.009306276, 0.009401238, 0.0094962,
      0.009591162, 0.009686124, 0.009781086, 0.009876048, 0.00997101,
      0.009113445, 0.009209376, 0.009305307, 0.009401238, 0.009497169,
      0.0095931, 0.009689031, 0.009784962, 0.009880893, 0.009976824,
      0.010072755, 0.0092055, 0.0093024, 0.0093993, 0.0094962, 0.0095931,
      0.00969, 0.0097869, 0.0098838, 0.0099807, 0.0100776, 0.0101745,
      0.009297555, 0.009395424, 0.009493293, 0.009591162, 0.009689031,
      0.0097869, 0.009884769, 0.009982638, 0.010080507, 0.010178376,
      0.010276245, 0.00938961, 0.009488448, 0.009587286, 0.009686124,
      0.009784962, 0.0098838, 0.009982638, 0.010081476, 0.010180314,
      0.010279152, 0.01037799, 0.009481665, 0.009581472, 0.009681279,
      0.009781086, 0.009880893, 0.0099807, 0.010080507, 0.010180314,
      0.010280121, 0.010379928, 0.010479735, 0.00957372, 0.009674496,
      0.009775272, 0.009876048, 0.009976824, 0.0100776, 0.010178376,
      0.010279152, 0.010379928, 0.010480704, 0.01058148, 0.009665775,
      0.00976752, 0.009869265, 0.00997101, 0.010072755, 0.0101745,
      0.010276245, 0.01037799, 0.010479735, 0.01058148, 0.010683225
      ), .Dim = c(11L, 11L))

x <- seq(from = 0.95, to = 1.05, by = 0.01)
y <- seq(from = 0.95, to = 1.05, by = 0.01)
run_f <- function(m, N) {
    m@gdata["upsilon"] * m@gdata["beta_t1"] * N
}
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
z_obs <- run_outer(x, y, model, upsilon ~ beta_t1, run_f, N = 3)

stopifnot(all(abs(z_obs - z_exp) < tol))

## Check missing 'x'
res <- tools::assertError(
    run_outer(y = y, model = model, formula = alpha ~ upsilon, FUN = run_f, N = 3))
check_error(res, "Missing 'x' argument.")

## Check non-numeric 'x'
res <- tools::assertError(
    run_outer(x = "x", y = y, model = model, formula = alpha ~ upsilon, FUN = run_f, N = 3))
check_error(res, "'x' argument is not numeric.")

## Check missing 'y'
res <- tools::assertError(
    run_outer(x = x, model = model, formula = alpha ~ upsilon, FUN = run_f, N = 3))
check_error(res, "Missing 'y' argument.")

## Check non-numeric 'y'
res <- tools::assertError(
    run_outer(x = x, y = "y", model = model, formula = alpha ~ upsilon, FUN = run_f, N = 3))
check_error(res, "'y' argument is not numeric.")

## Check missing 'model'
res <- tools::assertError(
    run_outer(x = x, y = y, formula = alpha ~ upsilon, FUN = run_f, N = 3))
check_error(res, "Missing 'model' argument.")

## Check non-SimInf_model 'model'
res <- tools::assertError(
    run_outer(x = x, y = y, model = "model", formula = alpha ~ upsilon, FUN = run_f, N = 3))
check_error(res, "'model' argument is not a 'SimInf_model'.")
