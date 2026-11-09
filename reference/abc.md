# Approximate Bayesian computation

Fit a model to data using the approximate Bayesian computation
sequential Monte Carlo (ABC SMC) algorithm of Toni and others (2009).
Each generation of the algorithm propagates a population of candidate
parameter sets (particles), drawn from the priors in the first
generation and from perturbed particles of the previous generation
thereafter, with the latter weighted by the ratio of their prior to
proposal densities. A particle is simulated by applying it to the model
and generating a stochastic trajectory with
[`run`](http://stewid.github.io/SimInf/reference/run.md). The particle
is retained if the distance between the summary statistics of the
simulated trajectory and the observed `data` is less than or equal to
the tolerance for the generation. As generations progress and the
tolerances decrease, the population of particles converges towards an
approximation of the posterior distribution.

## Usage

``` r
abc(
  model,
  priors = NULL,
  n_particles = NULL,
  n_init = NULL,
  distance = NULL,
  tolerance = NULL,
  data = NULL,
  verbose = FALSE,
  post_gen = NULL,
  init_model = NULL
)

# S4 method for class 'SimInf_model'
abc(
  model,
  priors = NULL,
  n_particles = NULL,
  n_init = NULL,
  distance = NULL,
  tolerance = NULL,
  data = NULL,
  verbose = FALSE,
  post_gen = NULL,
  init_model = NULL
)
```

## Arguments

- model:

  The `SimInf_model` object to generate data from.

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

  An integer `(>1)` specifying the number of particles to approximate
  the posterior with.

- n_init:

  Specify a positive integer (\>`n_particles`) to adaptively select a
  sequence of tolerances using the algorithm of Simola and others
  (2021). The initial tolerance is adaptively selected by sampling
  `n_init` draws from the prior and then retain the `n_particles`
  particles with the smallest distances. Note there must be enough
  initial particles to satisfactorily explore the parameter space, see
  Simola and others (2021). If the `tolerance` parameter is specified,
  then `n_init` must be `NULL`.

- distance:

  A function for calculating the summary statistics for a simulated
  trajectory. For each particle, the function must determine the
  distance and return that information. The first argument, `result`,
  passed to the `distance` function is the result from a `run` of the
  model with one trajectory attached to it. The second argument,
  `generation`, to `distance` is an integer with the generation of the
  particle(s). Further arguments that can passed to the `distance`
  function comes from the `data` argument in the `abc` function.
  Depending on the underlying model structure, data for one or more
  particles have been generated in each call to `distance`. If the
  `model` only contains one node and all the parameters to fit are in
  `ldata`, then that node will be replicated and each of the replicated
  nodes represent one particle in the trajectory (see ‘Examples’). On
  the other hand if the model contains multiple nodes or the parameters
  to fit are contained in `gdata`, then the trajectory in the `result`
  argument represents one particle. The function can return a numeric
  matrix (number of particles \\\times\\ number of summary statistics).
  Or, if the distance contains one summary statistic, a numeric vector
  with the length equal to the number of particles. Note that when using
  adaptive tolerance selection, only one summary statistic can be used,
  i.e., the function must return a matrix (number of particles
  \\\times\\ 1) or a numeric vector.

- tolerance:

  A numeric matrix (number of summary statistics \\\times\\ number of
  generations) where each column contains the tolerances for a
  generation and each row contains a sequence of gradually decreasing
  tolerances. Can also be a numeric vector if there is only one summary
  statistic. The tolerance determines the number of generations of
  ABC-SMC to run. If the `n_init` parameter is specified, then
  `tolerance` must be `NULL`.

- data:

  Optional data to be passed to the `distance` function. Default is
  `NULL`.

- verbose:

  prints diagnostic messages when `TRUE`. Default is `FALSE`.

- post_gen:

  An optional function that, if non-NULL, is applied after each
  completed generation. The function must accept one argument of type
  `SimInf_abc` with the current state of the fitting process. This
  function can be useful to, for example, save and inspect intermediate
  results.

- init_model:

  An optional function that, if non-NULL, is applied before running each
  proposal. The function must accept one argument of type `SimInf_model`
  with the current model of the fitting process. This function can be
  useful to specify the initial state of `u0` or `v0` of the model
  before running a trajectory with proposed parameters.

## Value

A `SimInf_abc` object.

## Details

The tolerances can either be supplied explicitly through the `tolerance`
argument, or determined adaptively using the algorithm of Simola and
others (2021), see `n_init`. The returned
[`SimInf_abc`](http://stewid.github.io/SimInf/reference/SimInf_abc-class.md)
object holds the accepted particles for each generation together with
their weights.

## References

T. Toni, D. Welch, N. Strelkowa, A. Ipsen, and M. P. H. Stumpf.
Approximate Bayesian computation scheme for parameter inference and
model selection in dynamical systems. *Journal of the Royal Society
Interface* **6**, 187–202, 2009.
[doi:10.1098/rsif.2008.0172](https://doi.org/10.1098/rsif.2008.0172)

U. Simola, J. Cisewski-Kehe, M. U. Gutmann, J. Corander. Adaptive
Approximate Bayesian Computation Tolerance Selection. *Bayesian
Analysis*, **16**(2), 397–423, 2021.
[doi:10.1214/20-BA1211](https://doi.org/10.1214/20-BA1211)

## See also

[`plot`](http://stewid.github.io/SimInf/reference/plot-SimInf_abc-method.md),
[`summary`](http://stewid.github.io/SimInf/reference/summary-SimInf_abc-method.md)
and
[`as.data.frame`](http://stewid.github.io/SimInf/reference/as.data.frame.SimInf_abc.md)
to inspect the fitted object,
[`continue_abc`](http://stewid.github.io/SimInf/reference/continue_abc.md)
to refine the fit with additional generations, and
[`run`](http://stewid.github.io/SimInf/reference/run.md) to generate
trajectories from the posterior.

[`continue_abc`](http://stewid.github.io/SimInf/reference/continue_abc.md).

## Examples

``` r
if (FALSE) { # \dontrun{
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
    ## trajectory count proportionally more against a proposal. Each
    ## proposal will be accepted or rejected depending on whether the
    ## distance is at most the tolerance for the current generation.
    residuals <- matrix(
        sim$I - data,
        nrow = length(data),
        ncol = n_nodes(result)
    )

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
} # }
```
