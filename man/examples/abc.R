## Let us consider an SIR model in a closed population with N = 100
## individuals of whom one is initially infectious and the rest are
## susceptible. First, generate one realisation (with a specified
## seed) from the model with known parameters \code{beta = 0.16} and
## \code{gamma = 0.077}. Then, use \code{abc} to infer the (known)
## parameters from the simulated data.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the SIR model and plot the number of infectious.
set.seed(22)
infectious <- trajectory(run(model), "I")$I
plot(infectious, type = "s")

## The distance function to accept or reject a proposal. Each node
## in the simulated trajectory (contained in the 'result' object)
## represents one proposal. The 'generation' argument is the current
## generation of proposals.
acceptFun <- function(result, generation, tol, ptol, ...) {
    ## Determine the tolerance for this generation.
    tol <- tol * (ptol)^(generation - 1)

    ## Extract the time-series of infectious in each node as a
    ## data.frame.
    sim <- trajectory(result, "I")

    ## Split the 'sim' data.frame by node and calculate the sum of the
    ## squared distance at each time-point for each node.
    dist <- tapply(sim$I, sim$node, function(sim_infectious) {
        sum((infectious - sim_infectious)^2)
    })

    ## Return TRUE or FALSE for each node depending on if the distance
    ## is less than the tolerance.
    abc_accept(dist < tol, tol)
}

## Fit the model parameters using ABC-SMC. The priors for the
## paramters are specified in the second argument using a formula
## notation. Here we use a uniform distribtion for each parameter with
## lower bound = 0 and upper bound = 1. Note that we use a low number
## particles here to keep the run-time of the example short.  In
## practice you would want to use many more to ensure better
## approximations.
fit <- abc(model = model,
           priors = c(beta ~ uniform(0, 1), gamma ~ uniform(0, 1)),
           ngen = 4,
           npart = 100,
           fn = acceptFun,
           tol = 5000,
           ptol = 0.9)

## Print a brief summary.
fit

## Display the ABC posterior distribution.
plot(fit)

## Run one more generation.
fit <- continue(fit, tol = 5000, ptol = 0.9)

plot(fit)
