% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimInf_model.R
\name{gdata}
\alias{gdata}
\alias{gdata,SimInf_model-method}
\title{Extract global data from a \code{SimInf_model} object}
\usage{
gdata(model)

\S4method{gdata}{SimInf_model}(model)
}
\arguments{
\item{model}{The \code{model} to get global data from.}
}
\value{
a numeric vector
}
\description{
The global data is a numeric vector that is common to all nodes.
The global data vector is passed as an argument to the transition
rate functions and the post time step function.
}
\examples{
## Create an SIR model
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:5, beta = 0.16, gamma = 0.077)

## Set 'beta' to a new value
gdata(model, "beta") <- 2

## Extract the global data vector that is common to all nodes
gdata(model)
}
