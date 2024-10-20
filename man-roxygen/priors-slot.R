##' @slot priors A \code{data.frame} containing the four columns
##'     \code{parameter}, \code{distribution}, \code{p1} and
##'     \code{p2}. The column \code{parameter} gives the name of the
##'     parameter referred to in the model. The column
##'     \code{distribution} contains the name of the prior
##'     distribution. Valid distributions are 'gamma', 'normal' or
##'     'uniform'. The column \code{p1} is a numeric vector with the
##'     first hyperparameter for each prior: 'gamma') shape,
##'     'lognormal') logmean, 'normal') mean, and 'uniform') lower
##'     bound. The column \code{p2} is a numeric vector with the
##'     second hyperparameter for each prior: 'gamma') rate,
##'     'lognormal') standard deviation on the log scale, 'normal')
##'     standard deviation, and 'uniform') upper bound.
