##' @slot U The result matrix with the number of individuals in each
##'     compartment in every node. \code{U[, j]} contains the number
##'     of individuals in each compartment at
##'     \code{tspan[j]}. \code{U[1:Nc, j]} contains the number of
##'     individuals in node 1 at \code{tspan[j]}. \code{U[(Nc + 1):(2
##'     * Nc), j]} contains the number of individuals in node 2 at
##'     \code{tspan[j]} etc. Integer matrix (\eqn{N_n N_c \times}
##'     \code{length(tspan)}).
