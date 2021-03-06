% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package_skeleton.R
\name{package_skeleton}
\alias{package_skeleton}
\title{Create a package skeleton from a \code{SimInf_model}}
\usage{
package_skeleton(
  model,
  name = NULL,
  path = ".",
  author = NULL,
  email = NULL,
  maintainer = NULL,
  license = "GPL-3"
)
}
\arguments{
\item{model}{The \code{model} \code{\linkS4class{SimInf_model}}
object with your model to create the package skeleton from.}

\item{name}{Character string with the package name. It should
contain only (ASCII) letters, numbers and dot, have at least
two characters and start with a letter and not end in a dot.
The package name is also used for the class name of the model
and the directory name of the package.}

\item{path}{Path to put the package directory in. Default is '.'
i.e. the current directory.}

\item{author}{Author of the package.}

\item{email}{Email of the package maintainer.}

\item{maintainer}{Maintainer of the package.}

\item{license}{License of the package. Default is 'GPL-3'.}
}
\value{
invisible \code{NULL}.
}
\description{
Describe your model in a logical way in R, then \code{mparse}
creates a \code{\linkS4class{SimInf_model}} object with your model
definition that can be installed as an add-on R package.
}
\references{
Read the \emph{Writing R Extensions} manual for more
    details.

Once you have created a \emph{source} package you need to install
it: see the \emph{R Installation and Administration} manual,
\code{\link{INSTALL}} and \code{\link{install.packages}}.
}
