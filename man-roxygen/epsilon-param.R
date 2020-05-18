##' @param epsilon A numeric vector with the incubation rate from
##'     exposed to infected where each node can have a different
##'     epsilon value. The vector must have length 1 or
##'     \code{nrow(u0)}.  If the vector has length 1, but the model
##'     contains more nodes, the epsilon value is repeated in all
##'     nodes.
