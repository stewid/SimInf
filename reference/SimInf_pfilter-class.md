# Class `SimInf_pfilter`

Storage class for the result of a bootstrap particle filter analysis.
The class holds the model with a sampled trajectory attached, the number
of particles used, the estimated log-likelihood, and the effective
sample size (ESS) at each time-point.

## Slots

- `model`:

  A
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  object with one complete sampled trajectory attached. At the end of
  the filtering, a single particle is drawn from the ensemble at the
  final time-point, and its full state evolution across all time-points
  is reconstructed.

- `n_particles`:

  An integer with the number of particles that was used at each
  time-step.

- `loglik`:

  The estimated log-likelihood.

- `ess`:

  A numeric vector with the effective sample size (ESS) at each
  time-point. The effective sample size is computed as
  \$\$\left(\sum\_{i=1}^N\\(w\_{t}^{i})^2\right)^{-1},\$\$ where
  \\w\_{t}^{i}\\ is the normalized weight of particle \\i\\ at time
  \\t\\.

## See also

[`pfilter`](http://stewid.github.io/SimInf/reference/pfilter.md) for
running a bootstrap particle filter and creating objects of this class,
[`logLik`](http://stewid.github.io/SimInf/reference/logLik-SimInf_pfilter-method.md)
for extracting the log-likelihood,
[`trajectory`](http://stewid.github.io/SimInf/reference/trajectory-SimInf_pfilter-method.md)
for extracting the filtered trajectory,
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence-SimInf_pfilter-method.md)
for computing prevalence from the filtered trajectory,
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for the underlying model class.
