## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
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

stopifnot(identical(
    SimInf:::match_compartments(compartments = ~S + I,
                                ok_combine = FALSE,
                                ok_lhs = FALSE,
                                U = c("S", "I", "R"),
                                V = "phi"),
    list(lhs = NULL,
         rhs = list(
             U = structure(
                 1:2,
                 .Names = c("S", "I"),
                 available_compartments = c("S", "I", "R")),
             V = structure(
                 integer(0),
                 .Names = character(0),
                 available_compartments = "phi")),
         condition = NULL)))

res <- assertError(
    SimInf:::match_compartments(compartments = ~S + I + phi,
                                ok_combine = FALSE,
                                ok_lhs = FALSE,
                                U = c("S", "I", "R"),
                                V = "phi"))
check_error(res, "Cannot combine data from different slots.")

stopifnot(identical(
    SimInf:::match_compartments(compartments = ~S + I + phi,
                                ok_combine = TRUE,
                                ok_lhs = FALSE,
                                U = c("S", "I", "R"),
                                V = "phi"),
    list(lhs = NULL,
         rhs = list(
             U = structure(
                 1:2,
                 .Names = c("S", "I"),
                 available_compartments = c("S", "I", "R")),
             V = structure(
                 c(phi = 1L),
                 available_compartments = "phi")),
         condition = NULL)))

res <- assertError(
    SimInf:::match_compartments(compartments = ~S + I + phi + test,
                                ok_combine = TRUE,
                                ok_lhs = FALSE,
                                U = c("S", "I", "R"),
                                V = "phi"))
check_error(res, "Non-existing compartment(s) in model: 'test'.")

stopifnot(identical(
    SimInf:::match_compartments(compartments = I ~ S + I + R,
                                ok_combine = FALSE,
                                ok_lhs = TRUE,
                                U = c("S", "I", "R"),
                                V = "phi"),
    list(lhs = list(
             U = structure(
                 c(I = 2L),
                 available_compartments = c("S", "I", "R")),
             V = structure(
                 integer(0),
                 .Names = character(0),
                 available_compartments = "phi")),
         rhs = list(
             U = structure(
                 1:3,
                 .Names = c("S", "I", "R"),
                 available_compartments = c("S", "I", "R")),
             V = structure(
                 integer(0),
                 .Names = character(0),
                 available_compartments = "phi")),
         condition = NULL)))

stopifnot(identical(
    SimInf:::match_compartments(compartments = I ~ S + I | R == 0,
                                ok_combine = FALSE,
                                ok_lhs = TRUE,
                                U = c("S", "I", "R"),
                                V = "phi"),
    list(lhs = list(
             U = structure(
                 c(I = 2L),
                 available_compartments = c("S", "I", "R")),
             V = structure(
                 integer(0),
                 .Names = character(0),
                 available_compartments = "phi")),
         rhs = list(
             U = structure(
                 1:2,
                 .Names = c("S", "I"),
                 available_compartments = c("S", "I", "R")),
             V = structure(
                 integer(0),
                 .Names = character(0),
                 available_compartments = "phi")),
         condition = "R == 0")))
