## Fit a model to data from an influenza in a boarding school in
## England (Anonymous. 1978. Influenza in a boarding school.
## British Medical Journal 1:578.)
obs <- data.frame(time = 1:15,
                  R1 = c(0, 1, 6, 26, 73, 222, 293, 258,
                         236, 191, 124, 69, 26, 11, 4))

## The distance function to accept or reject a proposal. Each node
## in the simulated trajectory (contained in the 'result' object)
## represents one proposal. The 'generation' argument is the current
## generation of proposals.
acceptFun <- function(result, generation, tol, ptol, ...) {
    ## Determine the tolerance for this generation.
    tol <- tol * (ptol)^(generation - 1)

    ## Extract the time-series for R1 for each node as a
    ## data.frame.
    sim <- trajectory(result, "R1")

    ## Split the 'sim' data.frame by node and calculate the sum
    ## of the squared distance at each time-point for each node.
    dist <- tapply(sim$R1, sim$node, function(sim_R1) {
        sum((obs$R1 - sim_R1)^2)
    })

    ## Return TRUE or FALSE for each node depending on if the
    ## distance is less than the tolerance.
    abc_accept(dist < tol, tol)
}

## Transitions
transitions <- c("S -> beta*S*I/(S+I+R1+R2) -> I",
                 "I -> gamma1*I -> R1",
                 "R1 -> gamma2*R1 -> R2")

## Specify the compartments.
compartments <- c("S", "I", "R1", "R2")

## Create the model.
model <- mparse(transitions = transitions,
                compartments = compartments,
                ldata = data.frame(beta = 0.1, gamma1 = 0.1, gamma2 = 0.1),
                u0 = data.frame(S = 762, I = 1, R1 = 0, R2 = 0),
                tspan = 1:15)

## Fit the model parameters using ABC-SMC. The priors for the paramters
## are specified in the second argument using a formula notation. Here
## we use a uniform distribtion for each parameter with lower bound = 0
## and upper bound = 5.
fit <- abc(model = model,
           priors = c(beta~U(0, 5), gamma1~U(0, 5), gamma2~U(0, 5)),
           ngen = 2,
           npart = 50,
           fn = acceptFun,
           tol = 100000,
           ptol = 0.5)

plot(fit)

fit <- continue(fit, tol = 100000, ptol = 0.5)

plot(fit)
