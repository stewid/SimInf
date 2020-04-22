##' @examples
##' \dontrun{
##' ## Use the model parser to create a 'SimInf_model' object that
##' ## expresses an SIR model, where 'beta' is the transmission rate
##' ## and 'gamma' is the recovery rate.
##' model  <- mparse(transitions = c("S -> beta*S*I/(S+I+R) -> I",
##'                                  "I -> gamma*I -> R"),
##'                  compartments = c("S", "I", "R"),
##'                  gdata = c(beta = 0.16, gamma = 0.077),
##'                  u0 = data.frame(S = 100, I = 1, R = 0),
##'                  tspan = 1:100)
##'
##' ## Run and plot the result
##' set.seed(22)
##' result <- run(model)
##' plot(result)
##' }
