##' @param U The result matrix with the number of individuals in each
##'     disease state in every node (\eqn{N_n N_c \times}
##'     \code{length(tspan)}).  \code{U[, j]} contains the number of
##'     individuals in each disease state at
##'     \code{tspan[j]}. \code{U[1:Nc, j]} contains the state of node
##'     \code{1} at \code{tspan[j]}. \code{U[(Nc + 1):(2 * Nc), j]}
##'     contains the state of node \code{2} at \code{tspan[j]} etc.
