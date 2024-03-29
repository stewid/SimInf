% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/u0.R
\name{u0<-}
\alias{u0<-}
\alias{u0<-,SimInf_model-method}
\title{Update the initial compartment state u0 in each node}
\usage{
u0(model) <- value

\S4method{u0}{SimInf_model}(model) <- value
}
\arguments{
\item{model}{The model to update the initial compartment state
\code{u0}.}

\item{value}{A \code{data.frame} with the initial state in each
node. Each row is one node, and the number of rows in
\code{u0} must match the number of nodes in \code{model}. Only
the columns in \code{u0} with a name that matches a
compartment in the \code{model} will be used.}
}
\description{
Update the initial compartment state u0 in each node
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

## Update u0 and run the model again
u0(model) <- data.frame(S = 990, I = 10, R = 0)
result <- run(model)
plot(result)
}
