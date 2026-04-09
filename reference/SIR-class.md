# Class SIR

Class to handle the SIR model. This class inherits from
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md),
meaning that SIR objects are fully compatible with all generic functions
defined for `SimInf_model`, such as
[`run`](http://stewid.github.io/SimInf/reference/run.md),
[`plot`](https://rdrr.io/r/graphics/plot.default.html),
[`trajectory`](http://stewid.github.io/SimInf/reference/trajectory.md),
and
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence.md).

## Details

The SIR model is a compartmental model for infectious diseases that
divides the population into three states: **S**usceptible, **I**nfected,
and **R**ecovered. It assumes that individuals gain permanent immunity
after recovery.

The model is defined by two state transitions: \$\$S \stackrel{\beta S I
/ N}{\longrightarrow} I\$\$ \$\$I \stackrel{\gamma I}{\longrightarrow}
R\$\$

where \\\beta\\ is the transmission rate, \\\gamma\\ is the recovery
rate, and \\N = S + I + R\\ is the total population size in each node.
Here, \\S\\, \\I\\, and \\R\\ represent the number of susceptible,
infected, and recovered individuals in that specific node.

## See also

[`SIR`](http://stewid.github.io/SimInf/reference/SIR.md) for creating an
SIR model object,
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for the parent class definition,
[`SEIR`](http://stewid.github.io/SimInf/reference/SEIR.md) for a model
including a latent period, and
[`SIS`](http://stewid.github.io/SimInf/reference/SIS.md) for a model
without immunity.
