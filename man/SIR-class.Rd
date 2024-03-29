% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIR.R
\docType{class}
\name{SIR-class}
\alias{SIR-class}
\title{Definition of the \acronym{SIR} model}
\description{
Class to handle the \acronym{SIR} \code{\link{SimInf_model}}.
}
\details{
The \acronym{SIR} model contains three compartments; number of
susceptible (S), number of infectious (I), and number of
recovered (R).  Moreover, it has two state transitions,
\deqn{S \stackrel{\beta S I / N}{\longrightarrow} I}{
  S -- beta S I / N --> I}
\deqn{I \stackrel{\gamma I}{\longrightarrow} R}{I -- gamma I --> R}
where \eqn{\beta} is the transmission rate, \eqn{\gamma} is the
recovery rate, and \eqn{N=S+I+R}.
}
\examples{
## Create an SIR model object.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the SIR model and plot the result.
set.seed(22)
result <- run(model)
plot(result)
}
