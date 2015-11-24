##' @section Beta:
##' The time dependent beta is divided into four intervals of the year
##' \preformatted{
##' where 0 <= day < 365 and
##' 0 <= end_t1 < end_t2 < end_t3 < end_t4 <= 365
##' or
##' 0 <= end_t4 < end_t1 < end_t2 < end_t3 < 364
##'
##' Case 1: end_t1 < end_t4
##' INTERVAL_1  INTERVAL_2       INTERVAL_3       INTERVAL_4       INTERVAL_1
##' [0, end_t1) [end_t1, end_t2) [end_t2, end_t3) [end_t3, end_t4) [end_t4, 365]
##'
##' Case 2: end_t4 < end_t1
##' INTERVAL_4  INTERVAL_1       INTERVAL_2       INTERVAL_3       INTERVAL_4
##' [0, end_t4) [end_t4, end_t1) [end_t1, end_t2) [end_t2, end_t3) [end_t3, 365)
##' }
