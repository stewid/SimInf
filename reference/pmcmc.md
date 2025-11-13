# Particle Markov chain Monte Carlo (PMCMC) algorithm

Particle Markov chain Monte Carlo (PMCMC) algorithm

## Usage

``` r
pmcmc(
  model,
  obs_process,
  data,
  priors,
  n_particles,
  n_iterations,
  theta = NULL,
  covmat = NULL,
  adaptmix = 0.05,
  adaptive = 100,
  post_proposal = NULL,
  init_model = NULL,
  post_particle = NULL,
  chain = NULL,
  verbose = getOption("verbose", FALSE)
)

# S4 method for class 'SimInf_model'
pmcmc(
  model,
  obs_process,
  data,
  priors,
  n_particles,
  n_iterations,
  theta = NULL,
  covmat = NULL,
  adaptmix = 0.05,
  adaptive = 100,
  post_proposal = NULL,
  init_model = NULL,
  post_particle = NULL,
  chain = NULL,
  verbose = getOption("verbose", FALSE)
)
```

## Arguments

- model:

  The model to simulate data from.

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

- priors:

  The priors for the parameters to fit. Each prior is specified with a
  formula notation, for example, `beta ~ uniform(0, 1)` specifies that
  beta is uniformly distributed between 0 and 1. Use
  [`c()`](https://rdrr.io/r/base/c.html) to provide more than one prior,
  for example, `c(beta ~ uniform(0, 1), gamma ~ normal(10, 1))`. The
  following distributions are supported: `gamma`, `lognormal`, `normal`
  and `uniform`. All parameters in `priors` must be only in either
  `gdata` or `ldata`.

- n_particles:

  An integer with the number of particles (\> 1) to use at each
  timestep.

- n_iterations:

  An integer specifying the number of iterations to run the PMCMC.

- theta:

  A named vector of initial values for the parameters of the model.
  Default is `NULL`, and then these are sampled from the prior
  distribution(s).

- covmat:

  A named numeric `(npars x npars)` matrix with covariances to use as
  initial proposal matrix. If left unspecified then defaults to
  `diag((theta/10)^2/npars)`.

- adaptmix:

  Mixing proportion for adaptive proposal. Must be a value between zero
  and one. Default is `adaptmix = 0.05`.

- adaptive:

  Controls when to start adaptive update. Must be greater or equal to
  zero. If `adaptive=0`, then adaptive update is not performed. Default
  is `adaptive = 100`.

- post_proposal:

  An optional function that, if non-`NULL`, is applied on the model
  after the proposal has been set for the model, but before running the
  particle filter. The function must accept one argument of type
  `SimInf_model` with the current model of the fitting process. This
  function can be useful to update, for example, `ldata` of the model
  before running a trajectory with proposed parameters. The function
  must return the model object which is then used in the particle
  filter.

- init_model:

  An optional function that, if non-NULL, is applied in the particle
  filter before running each proposal. The function must accept one
  argument of type `SimInf_model` with the current model of the fitting
  process. This function can be useful to specify the initial state of
  `u0` or `v0` of the model before running a trajectory with proposed
  parameters.

- post_particle:

  An optional function that, if non-NULL, is applied after each
  completed particle. The function must accept three arguments: 1) an
  object of `SimInf_pmcmc` with the current state of the fitting
  process, 2) an object `SimInf_pfilter` with the last particle and one
  filtered trajectory attached, and 3) an integer with the iteration in
  the fitting process. This function can be useful to, for example,
  monitor, save and inspect intermediate results. Note that the second
  `SimInf_pfilter` argument, is non-NULL only for the first particle in
  the chain, and for accepted particles.

- chain:

  An optional chain to start from. Must be a `data.frame` or an object
  that can be coerced to a `data.frame`. Only the columns in `chain`
  with a name that matches the names that will be used if this argument
  is not provided will be used. When this argument is provided,
  `n_iterations` can be 0. Additionally, when the `chain` argument is
  provided, then `theta` and `covmat` must be `NULL`.

- verbose:

  prints diagnostic messages when `TRUE`. The default is to retrieve the
  global option `verbose` and use `FALSE` if it is not set. When
  `verbose=TRUE`, information is printed every 100 iterations. For
  pmcmc, it is possible to get information every nth information by
  specifying `verbose=n`, for example, `verbose=1` or `verbose=10`.

## References

C. Andrieu, A. Doucet and R. Holenstein. Particle Markov chain Monte
Carlo methods. *Journal of the Royal Statistical Society, Series B*
**72**, 269–342, 2010.
[doi:10.1111/j.1467-9868.2009.00736.x](https://doi.org/10.1111/j.1467-9868.2009.00736.x)

G. O. Roberts and J. S. Rosenthal. Examples of adaptive MCMC. *Journal
of computational and graphical statistics*, **18**(2), 349–367, 2009.
[doi:10.1198/jcgs.2009.06134](https://doi.org/10.1198/jcgs.2009.06134)

## See also

[`continue_pmcmc`](http://stewid.github.io/SimInf/reference/continue_pmcmc.md).
