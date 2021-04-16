##' @param priors The priors for the parameters to fit. Each prior is
##'     specified with a formula notation, for example, \code{beta ~
##'     uniform(0, 1)} to specify that beta is uniformly distributed
##'     between 0 and 1. Use \code{c()} to provide more than one
##'     prior, for example, \code{c(beta ~ uniform(0, 1), gamma ~
##'     normal(10, 1)}. The following distributions are supported:
##'     \code{gamma}, \code{normal} and \code{uniform}.
