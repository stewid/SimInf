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
`mparse` serves as the primary interface for defining custom compartment
models, allowing users to describe transitions using a simple,
human-readable string syntax in R. The function then parses this
description, generates model-specific C code, and returns a
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object ready for simulation. This approach combines the flexibility of R
with the computational speed of compiled code, making it well-suited for
models with complex propensity functions, multiple compartments, or
node-specific parameters. Additionally, the package provides tools to
derive initial states directly from cleaned individual event logs.

After a model is created, a simulation is started with a call to the
[`run`](http://stewid.github.io/SimInf/reference/run.md) method and if
execution is successful, it returns a modified
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object with a single stochastic solution trajectory attached to it.

SimInf provides several utility functions to inspect simulated data, for
example, `show`, `summary` and `plot`. To facilitate custom analysis, it
provides the
[`trajectory,SimInf_model-method`](http://stewid.github.io/SimInf/reference/trajectory-SimInf_model-method.md)
and
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence.md)
methods.

One of our design goal was to make SimInf extendable and enable usage of
the numerical solvers from other R extension packages in order to
facilitate complex epidemiological research. To support this, SimInf has
functionality to generate the required C and R code from a model
specification, see
[`package_skeleton`](http://stewid.github.io/SimInf/reference/package_skeleton.md)

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
