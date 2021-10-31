##' @param obs_process Specification of the stochastic observation
##'     process. The \code{obs_process} can be specified as a
##'     \code{formula} if the model contains only one node and there
##'     is only one data point for each \code{time} in \code{data}.
##'     The left hand side of the formula must match a column name in
##'     the \code{data} data.frame and the right hand side of the
##'     formula is a character specifying the distribution of the
##'     observation process, for example, \code{Iobs ~ poisson(I)}.
##'     The following distributions are supported: \code{x ~
##'     binomial(size, prob)}, \code{x ~ poisson(rate)} and \code{x ~
##'     uniform(min, max)}. The observation process can also be a
##'     function to evaluate the probability density of the
##'     observations given the simulated states. The first argument
##'     passed to the \code{obs_process} function is the result from a
##'     run of the model and it contains one trajectory with simulated
##'     data for a time-point. The second argument to the
##'     \code{obs_process} function is a \code{data.frame} containing
##'     the rows for the specific time-point that the function is
##'     called for. Note that the function must return the log of the
##'     density.
