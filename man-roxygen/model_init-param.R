##' @param model_init An optional function that, if non-NULL, is
##'     applied before running each proposal. The function must accept
##'     one argument of type \code{SimInf_model} with the current
##'     model of the fitting process. This function can be useful to
##'     specify the initial state of u0 or v0 of the model before
##'     running a trajectory with proposed parameters.
