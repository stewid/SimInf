# Definition of the SIS model

Class to handle the SIS model. This class inherits from
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md),
meaning that SIS objects are fully compatible with all generic functions
defined for `SimInf_model`, such as
[`run`](http://stewid.github.io/SimInf/reference/run.md),
[`plot`](https://rdrr.io/r/graphics/plot.default.html),
[`trajectory`](http://stewid.github.io/SimInf/reference/trajectory.md),
and
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence.md).

## Details

The SIS model is a commonly used compartmental model for infectious
diseases where individuals do not gain permanent immunity after
recovery. Instead, they return to the susceptible state. It divides the
population into two states: **S**usceptible and **I**nfected.

The model is defined by two state transitions: \$\$S \stackrel{\beta S I
/ N}{\longrightarrow} I\$\$ \$\$I \stackrel{\gamma I}{\longrightarrow}
S\$\$

where \\\beta\\ is the transmission rate, \\\gamma\\ is the recovery
rate, and \\N = S + I\\ is the total population size in each node. Here,
\\S\\ and \\I\\ represent the number of susceptible and infected
individuals in that specific node.

## See also

[`SIS`](http://stewid.github.io/SimInf/reference/SIS.md) for creating an
SIS model object,
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for the parent class definition,
[`SIR`](http://stewid.github.io/SimInf/reference/SIR.md) for a model
with permanent immunity, and
[`SEIR`](http://stewid.github.io/SimInf/reference/SEIR.md) for a model
including a latent period.
