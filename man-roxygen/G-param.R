##' @param G Dependency graph that indicates the transition rates that
##'     need to be updated after a given state transition has occured.
##'     A non-zero entry in element \code{G[i, i]} indicates that
##'     transition rate \code{i} needs to be recalculated if the state
##'     transition \code{j} occurs. Sparse matrix (\eqn{Nt \times Nt})
##'     of object class
##'     \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}.
