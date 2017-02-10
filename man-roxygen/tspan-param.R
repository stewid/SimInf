##' @param tspan A vector (length >= 2) of increasing time points
##'     where the state of each node is to be returned. Can be either
##'     an \code{integer} or a \code{Date} vector. A \code{Date}
##'     vector is coerced to a numeric vector as days, where
##'     \code{tspan[1]} becomes the day of the year of the first year
##'     of \code{tspan}. The dates are added as names to the numeric
##'     vector.
