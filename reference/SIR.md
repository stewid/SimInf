# Create an SIR model

Create an SIR model to be used by the simulation framework.

## Usage

``` r
SIR(u0, tspan, events = NULL, beta = NULL, gamma = NULL)
```

## Arguments

- u0:

  A `data.frame` with the initial state in each node, i.e., the number
  of individuals in each compartment in each node when the simulation
  starts (see ‘Details’). The parameter `u0` can also be an object that
  can be coerced to a `data.frame`, e.g., a named numeric vector will be
  coerced to a one row `data.frame`.

- tspan:

  A vector (length \>= 1) of increasing time points where the state of
  each node is to be returned. Can be either an `integer` or a `Date`
  vector.

  - If `integer`: Represents absolute time steps.

  - If `Date`: Coerced to a numeric vector representing the **day of the
    year** (1–366) relative to the first date in the vector. The
    original `Date` objects are preserved as names for the numeric
    vector, facilitating time-series plotting.

- events:

  a `data.frame` with the scheduled events, see
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model.md).

- beta:

  A numeric vector with the transmission rate from susceptible to
  infected. Each node can have a different beta value. The vector must
  have length 1 or `nrow(u0)`. If the vector has length 1 but the model
  contains more nodes, the beta value is repeated for all nodes.

- gamma:

  A numeric vector with the recovery rate from infected to recovered.
  Each node can have a different gamma value. The vector must have
  length 1 or `nrow(u0)`. If the vector has length 1 but the model
  contains more nodes, the gamma value is repeated for all nodes.

## Value

A
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model.md)
of class `SIR`

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

The argument `u0` must be a `data.frame` with one row for each node with
the following columns:

- S:

  The number of susceptible individuals in each node

- I:

  The number of infected individuals in each node

- R:

  The number of recovered individuals in each node

## See also

[`SIR`](http://stewid.github.io/SimInf/reference/SIR-class.md) for the
class definition.
[`SEIR`](http://stewid.github.io/SimInf/reference/SEIR.md),
[`SIS`](http://stewid.github.io/SimInf/reference/SIS.md),
[`SISe`](http://stewid.github.io/SimInf/reference/SISe.md),
[`SISe3`](http://stewid.github.io/SimInf/reference/SISe3.md) and
[`SISe_sp`](http://stewid.github.io/SimInf/reference/SISe_sp.md) for
other predefined models.
[`mparse`](http://stewid.github.io/SimInf/reference/mparse.md) for
creating custom models.
[`run`](http://stewid.github.io/SimInf/reference/run.md) for running the
simulation.
[`trajectory`](http://stewid.github.io/SimInf/reference/trajectory.md),
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence.md)
and
[`plot,SimInf_model-method`](http://stewid.github.io/SimInf/reference/plot.md)
for post-processing and visualization.

## Examples

``` r
## For reproducibility, set the seed.
set.seed(22)

## Create an SIR model object.
model <- SIR(
  u0 = data.frame(S = 99, I = 1, R = 0),
  tspan = 1:100,
  beta = 0.16,
  gamma = 0.077
)

## Run the SIR model and plot the result.
result <- run(model)
plot(result)
```
