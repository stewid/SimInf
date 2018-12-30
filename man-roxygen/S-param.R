##' @param S Each column corresponds to a transition, and execution of
##'     state transition \code{j} amounts to adding the \code{S[, j]}
##'     to the state vector of the node where the state transition
##'     occurred.  Sparse matrix (\eqn{Nc \times Nt}) of object class
##'     \code{\linkS4class{dgCMatrix}}.
