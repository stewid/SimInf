##' @param post_proposal An optional function that, if
##'     non-\code{NULL}, is applied on the model after the proposal
##'     has been set for the model, but before running the particle
##'     filter. The function must accept one argument of type
##'     \code{SimInf_model} with the current model of the fitting
##'     process. This function can be useful to update, for example,
##'     \code{ldata} of the model before running a trajectory with
##'     proposed parameters. The function must return the model object
##'     which is then used in the particle filter.
