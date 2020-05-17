##' @param gamma A numeric vector with the recovery rate from infected
##'     to recovered where each node can have a different gamma
##'     value. The vector must have length 1 or \code{nrow(u0)}. If
##'     the vector has length 1, but the model contains more nodes,
##'     the beta value is repeated in all nodes.
