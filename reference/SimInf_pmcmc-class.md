# Class `SimInf_pmcmc`

Class `SimInf_pmcmc`

## Slots

- `model`:

  The `SimInf_model` object to estimate parameters in.

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

  Character vector (`"gdata"` or `"ldata"`) that determines if the
  ABC-SMC method estimates parameters in `model@gdata` (global data) or
  in `model@ldata` (local data).

- `pars`:

  An integer vector with the indices of the parameters in `target` that
  are being estimated.

- `n_particles`:

  An integer with the number of particles (\> 1) to use in the bootstrap
  particle filter.

- `data`:

  A `data.frame` holding the time series data for the observation
  process.

- `chain`:

  A matrix where each row contains `logPost`, `logLik`, `logPrior`,
  `accept`, and the `parameters` for each iteration.

- `covmat`:

  A named numeric `(npars x npars)` matrix with covariances to use as
  initial proposal matrix.

- `adaptmix`:

  Mixing proportion for adaptive proposal.

- `adaptive`:

  Controls when to start adaptive update.

## See also

[`pmcmc`](http://stewid.github.io/SimInf/reference/pmcmc.md) and
[`continue_pmcmc`](http://stewid.github.io/SimInf/reference/continue_pmcmc.md).
