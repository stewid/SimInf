# Definition of the SEIR model

Class to handle the SEIR model. This class inherits from
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md),
meaning that SEIR objects are fully compatible with all generic
functions defined for `SimInf_model`, such as
[`run`](http://stewid.github.io/SimInf/reference/run.md),
[`plot`](https://rdrr.io/r/graphics/plot.default.html),
[`trajectory`](http://stewid.github.io/SimInf/reference/trajectory.md),
and
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence.md).

## Details

The SEIR model extends the standard SIR model by adding an **E**xposed
(E) compartment for individuals who have been infected but are not yet
infectious. This accounts for the latent period of the disease.

The model is defined by three state transitions: \$\$S \stackrel{\beta S
I / N}{\longrightarrow} E\$\$ \$\$E \stackrel{\epsilon
E}{\longrightarrow} I\$\$ \$\$I \stackrel{\gamma I}{\longrightarrow}
R\$\$

where \\\beta\\ is the transmission rate, \\\epsilon\\ is the incubation
rate (inverse of the latent period), \\\gamma\\ is the recovery rate,
and \\N = S + E + I + R\\ is the total population size in each node.
Here, \\S\\, \\E\\, \\I\\, and \\R\\ represent the number of
susceptible, exposed, infected, and recovered individuals in that
specific node.

## See also

[`SEIR`](http://stewid.github.io/SimInf/reference/SEIR.md) for creating
an SEIR model object,
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for the parent class definition, and
[`SIR`](http://stewid.github.io/SimInf/reference/SIR.md) for the base
model without the latent period.
