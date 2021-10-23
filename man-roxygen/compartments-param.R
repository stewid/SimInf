##' @param compartments specify the names of the compartments to
##'     extract data from. The compartments can be specified as a
##'     character vector e.g. \code{compartments = c('S', 'I', 'R')},
##'     or as a formula e.g. \code{compartments = ~S+I+R} (see
##'     \sQuote{Examples}). Default (\code{compartments=NULL}) is to
##'     extract the number of individuals in each compartment i.e. the
##'     data from all discrete state compartments in the model. In
##'     models that also have continuous state variables e.g. the
##'     \code{SISe} model, they are also included.
