# Class `"SimInf_pmcmc"`

Class `"SimInf_pmcmc"`

## Slots

- `model`:

  The `SimInf_model` object to estimate parameters in.

- `priors`:

  A `data.frame` containing the four columns `parameter`,
  `distribution`, `p1` and `p2`. The column `parameter` gives the name
  of the parameter referred to in the model. The column `distribution`
  contains the name of the prior distribution. Valid distributions are
  'gamma', 'normal' or 'uniform'. The column `p1` is a numeric vector
  with the first hyperparameter for each prior: 'gamma') shape,
  'lognormal') logmean, 'normal') mean, and 'uniform') lower bound. The
  column `p2` is a numeric vector with the second hyperparameter for
  each prior: 'gamma') rate, 'lognormal') standard deviation on the log
  scale, 'normal') standard deviation, and 'uniform') upper bound.

- `target`:

  Character vector (`gdata` or `ldata`) that determines if the `pmcmc`
  method estimates parameters in `model@gdata` or in `model@ldata`.

- `pars`:

  Index to the parameters in `target`.

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
