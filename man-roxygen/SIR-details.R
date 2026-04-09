##' @details
##' The \acronym{SIR} model is a compartmental model for infectious
##' diseases that divides the population into three states:
##' \strong{S}usceptible, \strong{I}nfected, and \strong{R}ecovered.
##' It assumes that individuals gain permanent immunity after
##' recovery.
##'
##' The model is defined by two state transitions:
##' \deqn{S \stackrel{\beta S I / N}{\longrightarrow} I}{
##'   S -- beta S I / N --> I}
##' \deqn{I \stackrel{\gamma I}{\longrightarrow} R}{I -- gamma I --> R}
##'
##' where \eqn{\beta} is the transmission rate, \eqn{\gamma} is the
##' recovery rate, and \eqn{N = S + I + R} is the total population
##' size in each node. Here, \eqn{S}, \eqn{I}, and \eqn{R} represent
##' the number of susceptible, infected, and recovered individuals in
##' that specific node.
