##' @param formula A formula that specifies the compartments that
##'     define the cases with a disease or that have a specific
##'     characteristic (numerator), and the compartments that define
##'     the entire population of interest (denominator). The
##'     left-hand-side of the formula defines the cases, and the
##'     right-hand-side defines the population, for example,
##'     \code{I~S+I+R} in a \sQuote{SIR} model (see
##'     \sQuote{Examples}). The \code{.}  (dot) is expanded to all
##'     compartments, for example, \code{I~.}  is expanded to
##'     \code{I~S+I+R} in a \sQuote{SIR} model (see
##'     \sQuote{Examples}). The formula can also contain a condition
##'     (indicated by \code{|}) for each node and time step to further
##'     control the population to include in the calculation, for
##'     example, \code{I ~ . | R == 0} to calculate the prevalence
##'     when the recovered is zero in a \sQuote{SIR} model. The
##'     condition must evaluate to \code{TRUE} or \code{FALSE} in each
##'     node and time step. Note that if the denominator is zero, the
##'     prevalence is \code{NaN}.
