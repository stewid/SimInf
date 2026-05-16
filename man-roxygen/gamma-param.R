##' @param gamma A numeric vector with the recovery rate from infected
##'     to recovered. Each node can have a different gamma value.  The
##'     vector must have length 1 or \code{nrow(u0)}. If the vector
##'     has length 1 but the model contains more nodes, the gamma
##'     value is repeated for all nodes.
