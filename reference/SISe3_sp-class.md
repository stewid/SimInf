# Class SISe3_sp

Class to handle the SISe3_sp model. This class inherits from
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md),
meaning that SISe3_sp objects are fully compatible with all generic
functions defined for `SimInf_model`, such as
[`run`](http://stewid.github.io/SimInf/reference/run.md),
[`plot,SimInf_model-method`](http://stewid.github.io/SimInf/reference/plot.md),
[`trajectory`](http://stewid.github.io/SimInf/reference/trajectory.md),
and
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence.md).

## Details

The `SISe3_sp` model contains two compartments in three age categories:
**S**usceptible (\\S_1, S_2, S_3\\) and **I**nfected (\\I_1, I_2,
I_3\\). Additionally, it includes a continuous **environmental**
compartment (\\\varphi\\) to model shedding of a pathogen to the
environment. Moreover, it includes a spatial coupling of the
environmental contamination among proximal nodes to capture between-node
spread unrelated to moving infected individuals.

The model is defined by six state transitions:

\$\$S_1 \stackrel{\upsilon_1 \varphi S_1}{\longrightarrow} I_1\$\$
\$\$I_1 \stackrel{\gamma_1 I_1}{\longrightarrow} S_1\$\$ \$\$S_2
\stackrel{\upsilon_2 \varphi S_2}{\longrightarrow} I_2\$\$ \$\$I_2
\stackrel{\gamma_2 I_2}{\longrightarrow} S_2\$\$ \$\$S_3
\stackrel{\upsilon_3 \varphi S_3}{\longrightarrow} I_3\$\$ \$\$I_3
\stackrel{\gamma_3 I_3}{\longrightarrow} S_3\$\$

where the transition rate from susceptible to infected in age category
\\k\\ is proportional to the environmental contamination \\\varphi\\ and
the transmission rate \\\upsilon_k\\. The recovery rate \\\gamma_k\\
moves individuals from infected back to susceptible.

The environmental infectious pressure \\\varphi(t)\\ in each node
evolves according to:

\$\$\frac{d \varphi_i(t)}{dt} = \frac{\alpha \left(I\_{i,1}(t) +
I\_{i,2}(t) + I\_{i,3}(t)\right)}{N_i(t)} + \sum_k{\frac{\varphi_k(t)
N_k(t) - \varphi_i(t) N_i(t)}{N_i(t)} \cdot \frac{D}{d\_{ik}}} -
\beta(t) \varphi_i(t)\$\$

where \\\alpha\\ is the average shedding rate of the pathogen to the
environment per infected individual and \\N = S_1 + S_2 + S_3 + I_1 +
I_2 + I_3\\ the size of the node. Next comes the spatial coupling among
proximal nodes, where \\D\\ is the rate of the local spread and
\\d\_{ik}\\ the distance between holdings \\i\\ and \\k\\. The seasonal
decay and removal of the pathogen is captured by \\\beta(t)\\. The
environmental infectious pressure \\\varphi(t)\\ in each node is evolved
each time unit by the Euler forward method. The value of \\\varphi(t)\\
is saved at the time-points specified in `tspan`.

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

[`SISe3_sp`](http://stewid.github.io/SimInf/reference/SISe3_sp.md) for
creating an SISe3_sp model object and
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for the parent class definition.
