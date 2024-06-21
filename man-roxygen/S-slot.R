##' @slot S Each column corresponds to a state transition, and
##'     execution of state transition \code{j} amounts to adding the
##'     \code{S[, j]} column to the state vector \code{u[, i]} of node
##'     \emph{i} where the transition occurred. Sparse matrix (\eqn{Nc
##'     \times Nt}) of object class
##'     \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}.
