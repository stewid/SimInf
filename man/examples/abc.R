\dontrun{
## Let us consider an SIR model in a closed population with N = 100
## individuals of whom one is initially infectious and the rest are
## susceptible. First, generate one realisation (with a specified
## seed) from the model with known parameters beta=0.16 and
## gamma=0.077. Then, use abc to infer the (known) parameters from the
## simulated data.
model <- SIR(
  u0 = data.frame(
    S = 99,
    I = 1,
    R = 0),
  tspan = 1:100,
  beta = 0.16,
  gamma = 0.077
)

## Run the SIR model and plot the number of infectious.
set.seed(22)
infectious <- trajectory(run(model), "I")$I
plot(infectious, type = "s")

## The distance function to accept or reject a proposal. Each node in
## the simulated trajectory (contained in the 'result' object)
## represents one proposal.
distance <- function(result, generation, data) {
    ## Extract the time-series of infectious in each node as a
    ## data.frame, ordered so that time-points for the same node are
    ## contiguous.
    sim <- trajectory(result, "I")
    sim <- sim[order(sim$node, sim$time), ]

    ## One column per node: subtract the observed series row-wise. The
    ## residuals are squared and summed within each column, giving the
    ## sum of squared errors: larger deviations from the observed
    ## trajectory count proportionally more against a proposal.
    residuals <- matrix(
        sim$I - data,
        nrow = length(data),
        ncol = n_nodes(result)
    )

    ## Return the sum of squared errors for each node, i.e. one
    ## distance per proposal. Each proposal will be accepted or
    ## rejected depending on whether the distance is at most the
    ## tolerance for the current generation. The length of the
    ## returned vector must match the number of nodes in 'result'.
    colSums(residuals^2)
}

## Fit the model parameters using ABC-SMC and adaptive tolerance
## selection. The priors for the parameters are specified using a
## formula notation. Here we use a uniform distribution for each
## parameter with lower bound = 0 and upper bound = 1. Note that we
## use a low number of particles here to keep the run-time of the
## example short. In practice you would want to use many more to
## ensure better approximations.
fit <- abc(
  model = model,
  priors = c(
      beta ~ uniform(0, 1),
      gamma ~ uniform(0, 1)),
  n_particles = 100,
  n_init = 1000,
  distance = distance,
  data = infectious,
  verbose = TRUE
)

## Print a brief summary.
fit

## Display the ABC posterior distribution.
plot(fit)
}
