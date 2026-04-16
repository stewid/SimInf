##' @section Beta:
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
