##' @details
##' The \code{SISe} model contains two compartments:
##' \strong{S}usceptible (\eqn{S}) and \strong{I}nfected
##' (\eqn{I}). Additionally, it includes a continuous
##' \strong{environmental} compartment (\eqn{\varphi}) to model the
##' shedding of a pathogen to the environment.
##'
##' The model is defined by two state transitions:
##'
##' \deqn{S \stackrel{\upsilon \varphi S}{\longrightarrow} I}{ S --
##'   upsilon phi S --> I} \deqn{I \stackrel{\gamma
##'   I}{\longrightarrow} S}{ I -- gamma I --> S}
##'
##' where the transition rate from susceptible to infected is
##' proportional to the environmental contamination \eqn{\varphi} and
##' the transmission rate \eqn{\upsilon}. The recovery rate
##' \eqn{\gamma} moves individuals from infected back to susceptible.
##'
##' The environmental infectious pressure \eqn{\varphi(t)} in each
##' node evolves according to:
##'
##' \deqn{\frac{d\varphi(t)}{dt} = \frac{\alpha I(t)}{N(t)} - \beta(t)
##' \varphi(t) + \epsilon}{dphi/dt = alpha * sum(I) / N - beta(t) *
##' phi + epsilon}
##'
##' where:
##' \itemize{
##'   \item \eqn{\alpha} is the shedding rate per infected individual.
##'   \item \eqn{N(t) = S + I} is the total population size in the
##'     node.
##'   \item \eqn{\beta(t)} is the seasonal decay/removal rate, which
##'     varies throughout the year.
##'   \item \eqn{\epsilon} is the background infectious pressure.
##' }
##'
##' The environmental pressure is evolved using the Euler forward method
##' and saved at time points in \code{tspan}.
##'
##' \strong{Seasonal Decay (\eqn{\beta(t)}):}
##' The decay rate \eqn{\beta(t)} is piecewise constant, defined by four
##' intervals determined by the parameters \code{end_t1}, \code{end_t2},
##' \code{end_t3}, and \code{end_t4} (days of the year, where
##' \code{0 <= day < 365}). The year is divided into four intervals based
##' on the sorted order of these endpoints. The interval that wraps around
##' the year boundary (from the last endpoint to day 365, then from day 0
##' to the first endpoint) receives the same rate as the interval
##' preceding the first endpoint. Three orderings are supported:
##'
##' \strong{Case 1:} \code{end_t1 < end_t2 < end_t3 < end_t4}
##' \itemize{
##'   \item Interval 1: \code{[0, end_t1)} with rate \code{beta_t1}
##'   \item Interval 2: \code{[end_t1, end_t2)} with rate \code{beta_t2}
##'   \item Interval 3: \code{[end_t2, end_t3)} with rate \code{beta_t3}
##'   \item Interval 4: \code{[end_t3, end_t4)} with rate \code{beta_t4}
##'   \item Interval 1 (wrap-around): \code{[end_t4, 365)} with rate \code{beta_t1}
##' }
##'
##' \strong{Case 2:} \code{end_t3 < end_t4 < end_t1 < end_t2}
##' \itemize{
##'   \item Interval 3: \code{[0, end_t3)} with rate \code{beta_t3}
##'   \item Interval 4: \code{[end_t3, end_t4)} with rate \code{beta_t4}
##'   \item Interval 1: \code{[end_t4, end_t1)} with rate \code{beta_t1}
##'   \item Interval 2: \code{[end_t1, end_t2)} with rate \code{beta_t2}
##'   \item Interval 3 (wrap-around): \code{[end_t2, 365)} with rate \code{beta_t3}
##' }
##'
##' \strong{Case 3:} \code{end_t4 < end_t1 < end_t2 < end_t3}
##' \itemize{
##'   \item Interval 4: \code{[0, end_t4)} with rate \code{beta_t4}
##'   \item Interval 1: \code{[end_t4, end_t1)} with rate \code{beta_t1}
##'   \item Interval 2: \code{[end_t1, end_t2)} with rate \code{beta_t2}
##'   \item Interval 3: \code{[end_t2, end_t3)} with rate \code{beta_t3}
##'   \item Interval 4 (wrap-around): \code{[end_t3, 365)} with rate \code{beta_t4}
##' }
##'
##' These different orderings allow the model to handle seasonal patterns
##' where, for example, a winter peak crosses the year boundary.
