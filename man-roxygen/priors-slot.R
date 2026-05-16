##' @slot priors A \code{data.frame} defining the prior distributions
##'     for the parameters. It contains four columns:
##'     \itemize{
##'       \item \code{parameter}: The name of the parameter in the
##'       model.
##'       \item \code{distribution}: The prior distribution
##'         type. Valid values are \code{"gamma"}, \code{"lognormal"},
##'         \code{"normal"}, or \code{"uniform"}.
##'       \item \code{p1}: The first hyperparameter:
##'         \itemize{
##'           \item \code{"gamma"}: \emph{shape}
##'           \item \code{"lognormal"}: \emph{meanlog} (mean on the
##'           log scale)
##'           \item \code{"normal"}: \emph{mean}
##'           \item \code{"uniform"}: \emph{lower} bound
##'         }
##'       \item \code{p2}: The second hyperparameter:
##'         \itemize{
##'           \item \code{"gamma"}: \emph{rate}
##'           \item \code{"lognormal"}: \emph{sdlog} (standard
##'           deviation on the log scale)
##'           \item \code{"normal"}: \emph{sd} (standard deviation)
##'           \item \code{"uniform"}: \emph{upper} bound
##'         }
##'     }
