##' @slot U_sparse If the model was configured to write the solution
##'     to a sparse matrix
##'     (\code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}) the
##'     \code{U_sparse} contains the data and \code{U} is empty. The
##'     layout of the data in \code{U_sparse} is identical to
##'     \code{U}. Please note that \code{U_sparse} is numeric and
##'     \code{U} is integer.
