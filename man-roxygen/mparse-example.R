##' @examples
##' \dontrun{
##' ## Use the model parser to create a 'SimInf_mparse' object that
##' ## expresses a SIR model, where 'b' is the transmission rate and
##' ## 'g' is the recovery rate.
##' m <- mparse(transitions = c("S -> b*S*I/(S+I+R) -> I", "I -> g*I -> R"),
##'             compartments = c("S", "I", "R"),
##'             gdata = list(b = 0.16, g = 0.077))
##'
##' ## Initialize a 'SimInf_model' from the 'SimInf_mparse' object
##' u0 <- data.frame(S = 100, I = 1, R = 0)
##' model <- init(m, u0 = u0, tspan = 1:100)
##'
##' ## Run and plot the result
##' result <- run(model, threads = 1, seed = 22)
##' plot(result)
##' }
