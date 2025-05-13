## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 -- 2023 Stefan Widgren
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

## Define a tolerance
tol <- 1e-8

edges <- data.frame(
    from  = c(2, 3, 4, 1, 4, 5, 1, 3, 1, 3),
    to    = c(1, 1, 1, 2, 3, 3, 4, 4, 5, 5),
    rate  = c(0.2, 0.01, 0.79, 1, 0.2, 0.05, 0.2, 0.8, 0.2, 0.8),
    count = c(5, 5, 5, 50, 10, 10, 5, 5, 5, 5))

m_exp <- matrix(c(1, 0.2, 5, 2, 0.01, 5, 3, 0.79, 5, -1, 0, 1, 50, -1,
                  NaN, NaN, NaN, NaN, NaN, NaN, 3, 0.2, 10, 4, 0.05,
                  10, -1, NaN, NaN, NaN, 0, 0.2, 5, 2, 0.8, 5, -1,
                  NaN, NaN, NaN, 0, 0.2, 5, 2, 0.8, 5, -1, NaN, NaN,
                  NaN, -1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
                  NaN),
                nrow = 10L,
                ncol = 6L)

m_obs <- edge_properties_to_matrix(edges, 6)


stopifnot(identical(length(m_obs), length(m_exp)))
stopifnot(all(abs(m_obs[is.finite(m_obs)] - m_exp[is.finite(m_exp)]) < tol))

edges$from[1] <- NaN
res <- assertError(edge_properties_to_matrix(edges, 6))
check_error(res, "Values in 'edges' must be numeric and finite.")

edges$from[1] <- 2.2
res <- assertError(edge_properties_to_matrix(edges, 6))
check_error(res, "'edges' contain invalid 'from -> to' indices.")

edges$from[1] <- 0
res <- assertError(edge_properties_to_matrix(edges, 6))
check_error(res, "'edges' contain invalid 'from -> to' indices.")

edges$from[1] <- 2
edges <- rbind(edges[1, ], edges)
res <- assertError(edge_properties_to_matrix(edges, 6))
check_error(res, "'edges' contain duplicated properties.")
