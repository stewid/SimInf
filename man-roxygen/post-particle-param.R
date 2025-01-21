##' @param post_particle An optional function that, if non-NULL, is
##'     applied after each completed particle. The function must
##'     accept three arguments: 1) an object of \code{SimInf_pmcmc}
##'     with the current state of the fitting process, 2) an object
##'     \code{SimInf_pfilter} with the last particle and one filtered
##'     trajectory attached, and 3) an integer with the iteration in
##'     the fitting process. This function can be useful to, for
##'     example, monitor, save and inspect intermediate results. Note
##'     that the second \code{SimInf_pfilter} argument, is non-NULL
##'     only for the first particle in the chain, and for accepted
##'     particles.
