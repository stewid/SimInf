# A Framework for Data-Driven Stochastic Disease Spread Simulations

The SimInf package provides a flexible, high-performance framework for
data-driven spatio-temporal disease spread modeling. It is designed to
efficiently simulate disease transmission dynamics alongside population
demographics and dynamic contact networks.

## Details

The SimInf framework models infection dynamics within each subpopulation
(node) as continuous-time Markov chains (CTMC) using the Gillespie
stochastic simulation algorithm (SSA). Additionally, SimInf can
incorporate data—such as births, deaths, and movements—as scheduled
events. These events trigger at predefined time points and modify the
state of subpopulations by randomly sampling individuals from the
affected compartments. This capability allows simulations to be driven
by empirical records or synthetic scenarios while maintaining the
stochastic nature of the population dynamics.

The package supports both predefined models (e.g.,
[`SIR`](http://stewid.github.io/SimInf/reference/SIR.md),
[`SIS`](http://stewid.github.io/SimInf/reference/SIS.md)) and custom
model specifications via the
[`mparse`](http://stewid.github.io/SimInf/reference/mparse.md) function.
[`mparse`](http://stewid.github.io/SimInf/reference/mparse.md) serves as
the primary interface for defining custom compartment models, allowing
users to describe transitions using a simple, human-readable string
syntax in R. The function then parses this description, generates
model-specific C code, and returns a
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object ready for simulation. This approach combines the flexibility of R
with the computational speed of compiled code, making it well-suited for
models with complex propensity functions, multiple compartments, or
node-specific parameters. See the vignette "Getting started with mparse"
for a detailed tutorial on defining custom models.

To facilitate the use of real livestock data, SimInf provides
functionality for preparing individual event records (e.g., births,
deaths, and movements) for simulation. The
[`individual_events`](http://stewid.github.io/SimInf/reference/individual_events.md)
function cleans and validates raw livestock events, resolving
inconsistencies and structuring the data into the format required by the
simulation. This process transforms complex, individual-level logs into
the scheduled events that drive population demographics and movement
networks. See the vignette "Scheduled events" for a detailed guide on
processing livestock data and integrating it into disease spread models.

After a model is created, a simulation is executed using the
[`run`](http://stewid.github.io/SimInf/reference/run.md) method. Upon
successful completion,
[`run`](http://stewid.github.io/SimInf/reference/run.md) returns a new
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object containing the original configuration plus the simulated
stochastic trajectory.

To inspect and analyze the results, SimInf provides a suite of utility
functions:

- [`summary`](http://stewid.github.io/SimInf/reference/summary-SimInf_model-method.md)
  and
  [`show`](http://stewid.github.io/SimInf/reference/show-SimInf_model-method.md)
  for a quick overview of the model structure and simulation results.

- [`plot`](http://stewid.github.io/SimInf/reference/plot.md) for
  visualizing the time series of compartments and continuous state
  variables.

- [`trajectory`](http://stewid.github.io/SimInf/reference/trajectory-SimInf_model-method.md)
  for extracting the full time series data for custom analysis.

- [`prevalence`](http://stewid.github.io/SimInf/reference/prevalence.md)
  for calculating and summarizing disease prevalence across nodes and
  time.

These functions facilitate both rapid exploratory analysis and detailed
post-processing of simulation outcomes. See the vignette "Post-process
data in a trajectory" for a comprehensive tutorial on extracting and
analyzing simulation results.

Beyond simulation, the package provides functionality to fit models to
time series data using two Bayesian inference methods:

- Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC),
  implemented in
  [`abc`](http://stewid.github.io/SimInf/reference/abc.md), based on the
  approach by Toni and others (2009)
  [doi:10.1098/rsif.2008.0172](https://doi.org/10.1098/rsif.2008.0172) .

- Particle Markov Chain Monte Carlo (PMCMC), implemented in
  [`pmcmc`](http://stewid.github.io/SimInf/reference/pmcmc.md), based on
  the approach by Andrieu and others (2010)
  [doi:10.1111/j.1467-9868.2009.00736.x](https://doi.org/10.1111/j.1467-9868.2009.00736.x)
  .

Both methods enable parameter estimation in stochastic models where the
likelihood function is intractable, by using simulated data to estimate
the posterior distributions of model parameters.

Finally, one of the core design goals of SimInf is extensibility. The
package provides functionality to generate the required C and R code
from a model specification, enabling users to create custom R packages
that leverage the numerical solvers of SimInf. This facilitates complex
epidemiological research by allowing researchers to encapsulate their
specific model structures within a standard R package, ensuring
reproducibility and ease of sharing. See
[`package_skeleton`](http://stewid.github.io/SimInf/reference/package_skeleton.md)
for details on creating a new R package based on a custom model
specification.

## References

S. Widgren, P. Bauer, R. Eriksson and S. Engblom. SimInf: An R Package
for Data-Driven Stochastic Disease Spread Simulations. *Journal of
Statistical Software*, **91**(12), 1–42, 2019.
[doi:10.18637/jss.v091.i12](https://doi.org/10.18637/jss.v091.i12) . An
updated version of this paper is available as a vignette in the package.

## See also

Useful links:

- <https://github.com/stewid/SimInf>

- <https://stewid.github.io/SimInf/>

- Report bugs at <https://github.com/stewid/SimInf/issues>

## Author

**Maintainer**: Stefan Widgren <stefan.widgren@gmail.com>
([ORCID](https://orcid.org/0000-0001-5745-2284))

Authors:

- Robin Eriksson ([ORCID](https://orcid.org/0000-0002-4291-712X))

- Stefan Engblom ([ORCID](https://orcid.org/0000-0002-3614-1732))

- Pavol Bauer ([ORCID](https://orcid.org/0000-0003-4328-7171))

Other contributors:

- Thomas Rosendal ([ORCID](https://orcid.org/0000-0002-6576-9668))
  \[contributor\]

- Ivana Rodriguez Ewerlöf
  ([ORCID](https://orcid.org/0000-0002-9678-9813)) \[contributor\]

- Attractive Chaos (Author of 'kvec.h'.) \[copyright holder\]
