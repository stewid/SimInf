##' @param y Character vector or formula with the compartments in the
##'     model to include in the plot. Default includes all
##'     compartments in the model. Can also be a formula that
##'     specifies the compartments that define the cases with a
##'     disease or that have a specific characteristic (numerator),
##'     and the compartments that define the entire population of
##'     interest (denominator). The left-hand-side of the formula
##'     defines the cases, and the right-hand-side defines the
##'     population, for example, \code{I~S+I+R} in a \sQuote{SIR}
##'     model (see \sQuote{Examples}). The \code{.}  (dot) is expanded
##'     to all compartments, for example, \code{I~.}  is expanded to
##'     \code{I~S+I+R} in a \sQuote{SIR} model (see
##'     \sQuote{Examples}).
