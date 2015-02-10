## siminf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015  Stefan Engblom
## Copyright (C) 2015  Stefan Widgren
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

##' Generate a model for demonstration
##'
##' @param nodes Number of nodes in the model. Default is 1.
##' @param days Number of days to model. Default is 1000.
##' @param model The name of the model. Default is 'SISe'.
##' @return A model
##' @include external_events.r
##' @include SISe.r
##' @include SISe3.r
##' @export
demo_model <- function(nodes = 1,
                       days = 1000,
                       model = c("SISe"))
{
    ## Check 'nodes' argument
    if (!is.numeric(nodes))
        stop("'nodes' must be numeric.")
    if (!identical(length(nodes), 1L))
        stop("Length of 'nodes' must be one.")
    if (!is_wholenumber(nodes))
        stop("'nodes' must be integer")
    nodes <- as.integer(nodes)
    if (nodes[1] < 1)
        stop("'nodes' must be >= 1")

    ## Check 'days' argument
    if (!is.numeric(days))
        stop("'days' must be numeric.")
    if (!identical(length(days), 1L))
        stop("Length of 'days' must be one.")
    if (!is_wholenumber(days))
        stop("'days' must be integer")
    days <- as.integer(days)
    if (days[1] < 1)
        stop("'days' must be >= 1")

    ## Check 'model' argument
    model <- match.arg(model)

    if (identical(model, "SISe")) {
        model <- SISe(init    = data.frame(id = seq_len(nodes) - 1, S = 99, I = 1),
                      tspan   = seq_len(days) - 1,
                      events  = NULL,
                      phi     = rep(1, nodes),
                      upsilon = 0.017,
                      gamma   = 0.1,
                      alpha   = 1,
                      beta_q1 = -log(0.1) / 16,
                      beta_q2 = -log(0.1) / 16,
                      beta_q3 = -log(0.1) / 16,
                      beta_q4 = -log(0.1) / 16,
                      epsilon = 0)
    }

    model
}
