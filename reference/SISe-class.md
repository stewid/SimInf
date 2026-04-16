# Class SISe

Class to handle the SISe model. This class inherits from
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md),
meaning that SISe objects are fully compatible with all generic
functions defined for `SimInf_model`, such as
[`run`](http://stewid.github.io/SimInf/reference/run.md),
[`plot`](https://rdrr.io/r/graphics/plot.default.html),
[`trajectory`](http://stewid.github.io/SimInf/reference/trajectory.md),
and
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence.md).

## Details

The `SISe` model contains two compartments: **S**usceptible (\\S\\) and
**I**nfected (\\I\\). Additionally, it includes a continuous
**environmental** compartment (\\\varphi\\) to model the shedding of a
pathogen to the environment.

The model is defined by two state transitions:

\$\$S \stackrel{\upsilon \varphi S}{\longrightarrow} I\$\$ \$\$I
\stackrel{\gamma I}{\longrightarrow} S\$\$

where the transition rate from susceptible to infected is proportional
to the environmental contamination \\\varphi\\ and the transmission rate
\\\upsilon\\. The recovery rate \\\gamma\\ moves individuals from
infected back to susceptible.

The environmental infectious pressure \\\varphi(t)\\ in each node
evolves according to:

\$\$\frac{d\varphi(t)}{dt} = \frac{\alpha I(t)}{N(t)} - \beta(t)
\varphi(t) + \epsilon\$\$

where:

- \\\alpha\\ is the shedding rate per infected individual.

- \\N(t) = S + I\\ is the total population size in the node.

- \\\beta(t)\\ is the seasonal decay/removal rate, which varies
  throughout the year.

- \\\epsilon\\ is the background infectious pressure.

The environmental pressure is evolved using the Euler forward method and
saved at time points in `tspan`.

**Seasonal Decay (\\\beta(t)\\):** The decay rate \\\beta(t)\\ is
piecewise constant, defined by four intervals determined by the
parameters `end_t1`, `end_t2`, `end_t3`, and `end_t4` (days of the year,
where `0 <= day < 365`). The year is divided into four intervals based
on the sorted order of these endpoints. The interval that wraps around
the year boundary (from the last endpoint to day 365, then from day 0 to
the first endpoint) receives the same rate as the interval preceding the
first endpoint. Three orderings are supported:

**Case 1:** `end_t1 < end_t2 < end_t3 < end_t4`

- Interval 1: `[0, end_t1)` with rate `beta_t1`

- Interval 2: `[end_t1, end_t2)` with rate `beta_t2`

- Interval 3: `[end_t2, end_t3)` with rate `beta_t3`

- Interval 4: `[end_t3, end_t4)` with rate `beta_t4`

- Interval 1 (wrap-around): `[end_t4, 365)` with rate `beta_t1`

**Case 2:** `end_t3 < end_t4 < end_t1 < end_t2`

- Interval 3: `[0, end_t3)` with rate `beta_t3`

- Interval 4: `[end_t3, end_t4)` with rate `beta_t4`

- Interval 1: `[end_t4, end_t1)` with rate `beta_t1`

- Interval 2: `[end_t1, end_t2)` with rate `beta_t2`

- Interval 3 (wrap-around): `[end_t2, 365)` with rate `beta_t3`

**Case 3:** `end_t4 < end_t1 < end_t2 < end_t3`

- Interval 4: `[0, end_t4)` with rate `beta_t4`

- Interval 1: `[end_t4, end_t1)` with rate `beta_t1`

- Interval 2: `[end_t1, end_t2)` with rate `beta_t2`

- Interval 3: `[end_t2, end_t3)` with rate `beta_t3`

- Interval 4 (wrap-around): `[end_t3, 365)` with rate `beta_t4`

These different orderings allow the model to handle seasonal patterns
where, for example, a winter peak crosses the year boundary.

## See also

[`SISe`](http://stewid.github.io/SimInf/reference/SISe.md) for creating
an SISe model object and
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for the parent class definition.
