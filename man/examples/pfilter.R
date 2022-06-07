\dontrun{
## Let us consider an SIR model in a closed population with N = 100
## individuals of whom one is initially infectious and the rest are
## susceptible. First, generate one realisation (with a specified
## seed) from the model with known parameters \code{beta = 0.16} and
## \code{gamma = 0.077}.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = seq(1, 100, by = 3),
             beta = 0.16,
             gamma = 0.077)

## Run the SIR model to generate observed data for the number of
## infected individuals.
set.seed(22)
infected <- trajectory(run(model), "I")[, c("time", "I")]
colnames(infected) <- c("time", "Iobs")
}
