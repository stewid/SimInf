##' @param tspan A vector (length >= 1) of increasing time points
##'     where the state of each node is to be returned. Can be either
##'     an \code{integer} or a \code{Date} vector.
##'     \itemize{
##'       \item If \code{integer}: Represents the specific time points
##'         (e.g., days, hours) at which to record the state.
##'       \item If \code{Date}: Coerced to a numeric vector
##'         representing the \strong{day of the year} (1–366) relative
##'         to the first date in the vector. The original \code{Date}
##'         objects are preserved as names for the numeric vector,
##'         facilitating time-series plotting.
##'     }
