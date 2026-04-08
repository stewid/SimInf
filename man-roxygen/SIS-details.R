##' @details
##' The \acronym{SIS} model is a commonly used compartmental model for
##' infectious diseases where individuals do not gain permanent
##' immunity after recovery. Instead, they return to the susceptible
##' state.  It divides the population into two states:
##' \strong{S}usceptible and \strong{I}nfected.
##'
##' The model is defined by two state transitions:
##' \deqn{S \stackrel{\beta S I / N}{\longrightarrow} I}{S --> beta S
##'   I / N --> I}
##' \deqn{I \stackrel{\gamma I}{\longrightarrow} S}{I --> gamma I -->
##' S}
##'
##' where \eqn{\beta} is the transmission rate, \eqn{\gamma} is the
##' recovery rate, and \eqn{N = S + I} is the total population size in
##' each node. Here, \eqn{S} and \eqn{I} represent the number of
##' susceptible and infected individuals in that specific node.
