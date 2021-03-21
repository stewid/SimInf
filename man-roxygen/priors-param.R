##' @param priors The priors for the parameters to fit. Each prior is
##'     specified with a formula notation, for example, \code{beta ~
##'     U(0, 1)} to specify that beta is uniformly distributed between
##'     0 and 1. Use \code{c()} to provide more than one prior, for
##'     example, \code{c(beta ~ U(0, 1), gamma ~ N(10, 1)}. Gamma
##'     \code{G}, normal \code{N} and uniform \code{U} distributions
##'     are supported.
