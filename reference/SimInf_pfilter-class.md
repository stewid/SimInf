# Class `"SimInf_pfilter"`

Class `"SimInf_pfilter"`

## Slots

- `model`:

  A `SimInf_model` object with one filtered trajectory attached.

- `n_particles`:

  An integer with the number of particles that was used at each
  timestep.

- `loglik`:

  The estimated log likelihood.

- `ess`:

  A numeric vector with the effective sample size (ESS). The effective
  sample size is computed as
  \$\$\left(\sum\_{i=1}^N\\(w\_{t}^{i})^2\right)^{-1},\$\$ where
  \\w\_{t}^{i}\\ is the normalized weight of particle \\i\\ at time
  \\t\\.
