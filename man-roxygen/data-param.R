##' @param data A \code{data.frame} with the observed data. It must
##'     contain a \code{time} column with the observation time points
##'     of the data. The time points must be numeric and strictly
##'     increasing; rows that share the same time point are kept
##'     together. When \code{tspan} is expressed in days since a start
##'     date, the \code{time} column can also be provided as dates,
##'     which are converted to numeric time points. Additional columns
##'     hold the observed values and are made available, together with
##'     the time points, to the \code{obs_process}
##'     function. Observation times before the first time point in
##'     \code{tspan} are not allowed.
