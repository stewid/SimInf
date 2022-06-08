\dontrun{
## Let us consider an SIR model in a closed population with N = 100
## individuals of whom one is initially infectious and the rest are
## susceptible. First, generate one realisation (with a specified
## seed) from the model with known parameters 'beta = 0.16' and
## 'gamma = 0.077'. Then, use 'pfilter' to apply the bootstrap
## particle algorithm on the simulated data.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = seq(1, 100, by = 3),
             beta = 0.16,
             gamma = 0.077)

## Run the SIR model to generate simulated observed data for the
## number of infected individuals.
set.seed(22)
infected <- trajectory(run(model), "I")[, c("time", "I")]
colnames(infected) <- c("time", "Iobs")

## Use a Poison observation process for the infected individuals, such
## that 'Iobs ~ poison(I + 1e-6)'. A small constant '1e-6' is added to
## prevent numerical errors, since the simulated counts 'I' could be
## zero, which would result in the Poisson rate parameter being zero,
## which violates the conditions of the Poisson distribution. Use 1000
## particles.
pf <- pfilter(model,
              obs_process = Iobs ~ poisson(I + 1e-6),
              data = infected,
              npart = 1000)

## Print a brief summary.
pf

## Compare the number infected 'I' in the filtered trajectory with the
## infected 'Iobs' in the observed data.
plot(pf, ~I)
lines(Iobs ~ time, infected, col = "blue", lwd = 2, type = "s")
}
