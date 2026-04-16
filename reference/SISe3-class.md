# Class SISe3

Class to handle the SISe3 model. This class inherits from
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md),
meaning that SISe3 objects are fully compatible with all generic
functions defined for `SimInf_model`, such as
[`run`](http://stewid.github.io/SimInf/reference/run.md),
[`plot`](https://rdrr.io/r/graphics/plot.default.html),
[`trajectory`](http://stewid.github.io/SimInf/reference/trajectory.md),
and
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence.md).

## Details

The `SISe3` model contains two compartments in three age categories:
**S**usceptible (\\S_1, S_2, S_3\\) and **I**nfected (\\I_1, I_2,
I_3\\). Additionally, it includes a continuous **environmental**
compartment (\\\varphi\\) to model the shedding of a pathogen to the
environment.

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

\$\$\frac{d\varphi(t)}{dt} = \frac{\alpha \left(I_1(t) + I_2(t) +
I_3(t)\right)}{N(t)} - \beta(t) \varphi(t) + \epsilon\$\$

where:

- \\\alpha\\ is the shedding rate per infected individual.

- \\N(t) = S_1 + S_2 + S_3 + I_1 + I_2 + I_3\\ is the total population
  size in the node.

- \\\beta(t)\\ is the seasonal decay/removal rate, which varies
  throughout the year.

- \\\epsilon\\ is the background infectious pressure.

The environmental pressure is evolved using the Euler forward method and
saved at time points in `tspan`.
