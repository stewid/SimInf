##' @details
##'
##' The \acronym{SEIR} model extends the standard \acronym{SIR} model
##' by adding an \strong{E}xposed (E) compartment for individuals who
##' have been infected but are not yet infectious. This accounts for
##' the latent period of the disease.
##'
##' The model is defined by three state transitions:
##' \deqn{S \stackrel{\beta S I / N}{\longrightarrow} E}{ S -- beta S
##'   I / N --> E}
##' \deqn{E \stackrel{\epsilon E}{\longrightarrow} I}{E -- epsilon E
##' --> I}
##' \deqn{I \stackrel{\gamma I}{\longrightarrow} R}{I -- gamma I -->
##' R}
##'
##' where \eqn{\beta} is the transmission rate, \eqn{\epsilon} is the
##' incubation rate (inverse of the latent period), \eqn{\gamma} is
##' the recovery rate, and \eqn{N = S + E + I + R} is the total
##' population size in each node. Here, \eqn{S}, \eqn{E}, \eqn{I}, and
##' \eqn{R} represent the number of susceptible, exposed, infected,
##' and recovered individuals in that specific node.
