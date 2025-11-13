# Bootstrap particle filter

The bootstrap filtering algorithm. Systematic resampling is performed at
each observation.

## Usage

``` r
pfilter(model, obs_process, data, n_particles, init_model = NULL)

# S4 method for class 'SimInf_model'
pfilter(model, obs_process, data, n_particles, init_model = NULL)
```

## Arguments

- model:

  The `SimInf_model` object to simulate data from.

- obs_process:

  Specification of the stochastic observation process. The `obs_process`
  can be specified as a `formula` if the model contains only one node
  and there is only one data point for each `time` in `data`. The left
  hand side of the formula must match a column name in the `data`
  data.frame and the right hand side of the formula is a character
  specifying the distribution of the observation process, for example,
  `Iobs ~ poisson(I)`. The following distributions are supported:
  `x ~ binomial(size, prob)`, `x ~ poisson(rate)` and
  `x ~ uniform(min, max)`. The observation process can also be a
  function to evaluate the probability density of the observations given
  the simulated states. The first argument passed to the `obs_process`
  function is the result from a run of the model and it contains one
  trajectory with simulated data for a time-point, where the trajectory
  contains `n_particles` replicates, see
  [`trajectory,SimInf_model-method`](http://stewid.github.io/SimInf/reference/trajectory-SimInf_model-method.md).
  The second argument to the `obs_process` function is a `data.frame`
  containing the rows for the specific time-point that the function is
  called for. Note that the function must return the log of the density.

- data:

  A `data.frame` holding the time series data.

- n_particles:

  An integer with the number of particles (\> 1) to use at each
  timestep.

- init_model:

  An optional function that, if non-NULL, is applied before running each
  proposal. The function must accept one argument of type `SimInf_model`
  with the current model of the fitting process. This function can be
  useful to specify the initial state of `u0` or `v0` of the model
  before running a trajectory with proposed parameters.

## Value

A `SimInf_pfilter` object.

## References

N. J. Gordon, D. J. Salmond, and A. F. M. Smith. Novel Approach to
Nonlinear/Non-Gaussian Bayesian State Estimation. *Radar and Signal
Processing, IEE Proceedings F*, **140**(2) 107–113, 1993.
[doi:10.1049/ip-f-2.1993.0015](https://doi.org/10.1049/ip-f-2.1993.0015)

## Examples

``` r
if (FALSE) { # \dontrun{
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
              n_particles = 1000)

## Print a brief summary.
pf

## Compare the number infected 'I' in the filtered trajectory with the
## infected 'Iobs' in the observed data.
plot(pf, ~I)
lines(Iobs ~ time, infected, col = "blue", lwd = 2, type = "s")
} # }
```
