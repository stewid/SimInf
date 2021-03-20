##' @slot priors A \code{data.frame} containing the four columns
##'     \code{parameter}, \code{distribution}, \code{p1} and
##'     \code{p2}. The column \code{parameter} gives the name of the
##'     parameter referred to in the model. The column
##'     \code{distribution} contains a letter indicating the prior
##'     distribution. Valid letters are 'G' (gamma), 'N' (normal) or
##'     'U' (uniform). The column \code{p1} is a numeric vector with
##'     the first hyperparameter for each prior: 'G') shape, 'N')
##'     mean, and 'U') lower bound. The column \code{p2} is a numeric
##'     vector with the second hyperparameter for each prior: 'G')
##'     rate, 'N') standard deviation, and 'U') upper bound.
