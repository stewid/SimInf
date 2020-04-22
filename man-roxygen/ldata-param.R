##' @param ldata local data for the nodes. Can either be specified as
##'     a \code{data.frame} with one row per node. Or as a matrix
##'     where each column \code{ldata[, j]} contains the local data
##'     vector for the node \code{j}. The local data vector is passed
##'     as an argument to the transition rate functions and the post
##'     time step function.
