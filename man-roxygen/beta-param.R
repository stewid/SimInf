##' @param beta A numeric vector with the transmission rate from
##'     susceptible to infected where each node can have a different
##'     beta value. The vector must have length 1 or \code{nrow(u0)}.
##'     If the vector has length 1, but the model contains more nodes,
##'     the beta value is repeated in all nodes.
