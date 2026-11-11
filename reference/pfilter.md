# Bootstrap particle filter

Estimate the hidden trajectory of a stochastic model from observations
using the bootstrap particle filter (Gordon et al., 1993).

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
  function is the result from a run of the model for the current time
  interval: a trajectory with one replicate per particle for that
  time-point, see
  [`trajectory,SimInf_model-method`](http://stewid.github.io/SimInf/reference/trajectory.md).
  The second argument to the `obs_process` function is a `data.frame`
  containing the rows in `data` for the time-point that the function is
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

A `SimInf_pfilter` object with the filtered trajectory: a single
trajectory sampled from the particle system (traced backward through the
particle genealogy), together with the estimated log-likelihood of the
data and the effective sample size at each observation. See
[`SimInf_pfilter-class`](http://stewid.github.io/SimInf/reference/SimInf_pfilter-class.md)
for details.

## Details

The particle filter approximates the distribution of the latent state
conditional on the observed data. A large number of independent
trajectories ('particles') of the model are simulated forward in time.
At each observation, every particle is assigned a weight given by the
likelihood of the observation under that particle. Particles with a poor
fit to the data therefore receive a low weight. Particles are resampled
proportionally to their weights using systematic resampling, and the
resulting particles are propagated forward to the next observation. This
focuses the computational effort on trajectories that are consistent
with the data and yields a weighted approximation of the filtered state
distribution at each observation time. The accuracy of the approximation
improves as the number of particles increases. After the final
observation, a single trajectory is sampled from the particle system and
traced backward through the ancestral lineages. The effective sample
size is recorded at each observation, before resampling, as a measure of
how unevenly the weights were distributed at that step.

## References

N. J. Gordon, D. J. Salmond, and A. F. M. Smith. Novel Approach to
Nonlinear/Non-Gaussian Bayesian State Estimation. *Radar and Signal
Processing, IEE Proceedings F*, **140**(2) 107–113, 1993.
[doi:10.1049/ip-f-2.1993.0015](https://doi.org/10.1049/ip-f-2.1993.0015)

## Examples

``` r
## For reproducibility, specify a seed and the number of threads to
## use.
set.seed(22)
set_num_threads(1)

## Let us consider an SIR model in a closed population with 'N = 100'
## individuals of whom one is initially infectious and the rest
## are susceptible. First, generate one realisation from the model
## with known parameters 'beta = 0.16' and 'gamma = 0.077'. Specify
## 'tspan' to record the process every third day. Then, use
## 'pfilter' to apply the bootstrap particle filter algorithm on
## the simulated data.
model <- SIR(
  u0 = data.frame(
    S = 99,
    I = 1,
    R = 0),
  tspan = seq(1, 100, by = 3),
  beta = 0.16,
  gamma = 0.077
)

## Run the SIR model to generate the latent number of infected
## individuals, then simulate noisy observations from them using the
## Poisson observation model.
infected <- trajectory(run(model))
infected$Iobs <- rpois(nrow(infected), infected$I)

## Use a Poisson observation process for the number of infected
## individuals: Iobs ~ poisson(I + 1e-6). A small constant '1e-6' is
## added because the latent counts 'I' can be zero, in which case the
## Poisson density is zero for any positive observation. In the
## particle filter this would give the particle zero likelihood, so
## the small constant keeps the weights finite.

## Use a low number of particles here to keep the run-time of
## the example short. In practice, you may need to increase the
## number of particles for a smoother filtered trajectory.
pf <- pfilter(
  model,
  obs_process = Iobs ~ poisson(I + 1e-6),
  data = infected[, c("time", "Iobs")],
  n_particles = 100
)

## Print a brief summary.
pf
#> Number of particles: 100
#> Log-likelihood: -87.932284

## Compare the number of infected individuals 'I' in the filtered
## trajectory with the noisy observations 'Iobs' and the latent
## trajectory 'I'.
plot(pf, ~I, ylim = c(0, 30))

lines(
  Iobs ~ time,
  infected,
  col = "blue",
  lwd = 2,
  type = "s",
  lty = "dotted"
)

lines(
  I ~ time,
  infected,
  col = "red",
  lwd = 2,
  type = "s",
  lty = "dashed"
)

## Add a legend.
legend(
  "topright",
  legend = c("Filtered", "Observed", "Latent"),
  col = c("black", "blue", "red"),
  lty = c("solid", "dotted", "dashed"),
  lwd = 2
)
```
