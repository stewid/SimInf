##' @section Beta:
##' The time dependent beta is divided into four intervals of the year
##' \preformatted{
##' where 0 <= day < 365
##'
##' Case 1: END_1 < END_2 < END_3 < END_4
##' INTERVAL_1 INTERVAL_2     INTERVAL_3     INTERVAL_4     INTERVAL_1
##' [0, END_1) [END_1, END_2) [END_2, END_3) [END_3, END_4) [END_4, 365)
##'
##' Case 2: END_3 < END_4 < END_1 < END_2
##' INTERVAL_3 INTERVAL_4     INTERVAL_1     INTERVAL_2     INTERVAL_3
##' [0, END_3) [END_3, END_4) [END_4, END_1) [END_1, END_2) [END_2, 365)
##'
##' Case 3: END_4 < END_1 < END_2 < END_3
##' INTERVAL_4 INTERVAL_1     INTERVAL_2     INTERVAL_3     INTERVAL_4
##' [0, END_4) [END_4, END_1) [END_1, END_2) [END_2, END_3) [END_3, 365)
##' }
