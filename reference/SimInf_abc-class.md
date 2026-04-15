# Class `SimInf_abc`

Storage class for the results of an Approximate Bayesian Computation
(ABC) parameter estimation using Sequential Monte Carlo (SMC). The
`SimInf_abc` class holds the model definition, prior distributions,
accepted parameter values (particles), weights, distances, and
convergence diagnostics.

## Slots

- `model`:

  A
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  object containing the model structure (transitions, compartments,
  etc.) for which parameters are being estimated.

- `priors`:

  A `data.frame` defining the prior distributions for the parameters. It
  contains four columns:

  - `parameter`: The name of the parameter in the model.

  - `distribution`: The prior distribution type. Valid values are
    `"gamma"`, `"lognormal"`, `"normal"`, or `"uniform"`.

  - `p1`: The first hyperparameter:

    - `"gamma"`: *shape*

    - `"lognormal"`: *meanlog* (mean on the log scale)

    - `"normal"`: *mean*

    - `"uniform"`: *lower* bound

  - `p2`: The second hyperparameter:

    - `"gamma"`: *rate*

    - `"lognormal"`: *sdlog* (standard deviation on the log scale)

    - `"normal"`: *sd* (standard deviation)

    - `"uniform"`: *upper* bound

- `target`:

  Character vector (`gdata` or `ldata`) that determines if the ABC-SMC
  method estimates parameters in `model@gdata` or in `model@ldata`.

- `pars`:

  Index to the parameters in `target`.

- `nprop`:

  An integer vector with the number of simulated proposals in each
  generation.

- `fn`:

  A function for calculating the summary statistics for the simulated
  trajectory and determine the distance for each particle, see
  [`abc`](http://stewid.github.io/SimInf/reference/abc.md) for more
  details.

- `tolerance`:

  A numeric matrix (number of summary statistics \\\times\\ number of
  generations) where each column contains the tolerances for a
  generation and each row contains a sequence of gradually decreasing
  tolerances.

- `x`:

  A numeric array (number of particles \\\times\\ number of parameters
  \\\times\\ number of generations) with the parameter values for the
  accepted particles in each generation. Each row is one particle.

- `weight`:

  A numeric matrix (number of particles \\\times\\ number of
  generations) with the weights for the particles `x` in the
  corresponding generation.

- `distance`:

  A numeric array (number of particles \\\times\\ number of summary
  statistics \\\times\\ number of generations) with the distance for the
  particles `x` in each generation. Each row contains the distance for a
  particle and each column contains the distance for a summary
  statistic.

- `ess`:

  A numeric vector with the effective sample size (ESS) in each
  generation. The effective sample size is computed as
  \$\$\left(\sum\_{i=1}^N\\(w\_{g}^{(i)})^2\right)^{-1},\$\$ where
  \\w\_{g}^{(i)}\\ is the normalized weight of particle \\i\\ in
  generation \\g\\.

- `init_model`:

  An optional function that, if non-NULL, is applied before running each
  proposal. The function must accept one argument of type `SimInf_model`
  with the current model of the fitting process. This function can be
  useful to specify the initial state of `u0` or `v0` of the model
  before running a trajectory with proposed parameters.

## See also

[`abc`](http://stewid.github.io/SimInf/reference/abc.md) and
[`continue_abc`](http://stewid.github.io/SimInf/reference/continue_abc.md).
