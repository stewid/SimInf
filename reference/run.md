# Run the SimInf stochastic simulation algorithm

Simulate a stochastic trajectory from a SimInf model. The function
compiles and loads model-specific C code (if not already compiled),
initializes the solver, and advances the simulation in continuous time
from the first to the last time point in `tspan`, recording the state of
the system at each of these time points. Sparse output can be requested
with
[`punchcard<-`](http://stewid.github.io/SimInf/reference/punchcard-set.md)
to store only selected time points or compartments.

## Usage

``` r
run(model, ...)

# S4 method for class 'SimInf_model'
run(model, ...)

# S4 method for class 'SEIR'
run(model, ...)

# S4 method for class 'SIR'
run(model, ...)

# S4 method for class 'SIS'
run(model, ...)

# S4 method for class 'SISe'
run(model, ...)

# S4 method for class 'SISe3'
run(model, ...)

# S4 method for class 'SISe3_sp'
run(model, ...)

# S4 method for class 'SISe_sp'
run(model, ...)

# S4 method for class 'SimInf_abc'
run(model, ...)
```

## Arguments

- model:

  The
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  (or
  [`SimInf_abc`](http://stewid.github.io/SimInf/reference/SimInf_abc-class.md))
  object to run.

- ...:

  Optional arguments that affect the simulation:

  solver

  :   Character string specifying the numerical solver. Either `"ssm"`
      (default) or `"aem"`. If not provided, the `"ssm"` solver is used.
      See details.

  seed

  :   Numeric or integer specifying the random seed for the simulation.
      If not provided, a seed is randomly sampled from the current R RNG
      state.

  crn

  :   Logical. If `TRUE`, use common random numbers for the `"mssm"`
      solver by assigning each node and replicate a separate random
      number stream, enabling variance reduction across simulation runs.
      If not provided, `crn` is `FALSE`.

  rng_type

  :   Character string specifying the random number generator algorithm.
      Either `"mt19937"` or `"taus2"`. If not provided, `"mt19937"` is
      used.

## Value

The model object with a single stochastic trajectory attached.

## Details

The solver uses a split-step method: for each unit of time, it first
integrates the continuous-time Markov chain within each node using
direct SSA (Gillespie's algorithm), then processes scheduled events
(exit, enter, internal transfer, and external transfer), and finally
calls the post time step function to update continuous state variables;
see Widgren and others (2019) for details.

Two numerical solvers are available. The default solver is `"ssm"`
(split-step method), which becomes `"mssm"` (multi-model split-step
method) automatically when the model contains multiple replicates
(`replicates > 1`). The alternative solver is `"aem"` (all events
method), which assigns each reaction channel its own random number
stream for consistent operational times across simulations; see Bauer
and Engblom (2015). The `"aem"` solver cannot be used with replicated
models.

For a
[`SimInf_abc`](http://stewid.github.io/SimInf/reference/SimInf_abc-class.md)
object, `run()` samples one parameter set (a particle) from the last
generation of an approximate Bayesian computation sequential Monte Carlo
(ABC SMC) analysis and applies it to the underlying model before
simulating a stochastic trajectory. The particle is drawn uniformly at
random from the particles in the final generation, such that repeated
calls to `run()` yield simulations distributed according to the
approximated posterior distribution.

## References

S. Widgren, P. Bauer, R. Eriksson and S. Engblom. SimInf: An R Package
for Data-Driven Stochastic Disease Spread Simulations. *Journal of
Statistical Software*, **91**(12), 1–42, 2019.
[doi:10.18637/jss.v091.i12](https://doi.org/10.18637/jss.v091.i12) . An
updated version of this paper is available as a vignette in the package.

P. Bauer, S. Engblom and S. Widgren. Fast Event-Based Epidemiological
Simulations on National Scales. *International Journal of High
Performance Computing Applications*, **30**(4), 438–453, 2016. doi:
10.1177/1094342016635723

P. Bauer and S. Engblom. Sensitivity Estimation and Inverse Problems in
Spatial Stochastic Models of Chemical Kinetics. In: A. Abdulle, S.
Deparis, D. Kressner, F. Nobile and M. Picasso (eds.), *Numerical
Mathematics and Advanced Applications - ENUMATH 2013*, pp. 519–527,
Lecture Notes in Computational Science and Engineering, vol 103.
Springer, Cham, 2015.
[doi:10.1007/978-3-319-10705-9_51](https://doi.org/10.1007/978-3-319-10705-9_51)

## See also

[`trajectory`](http://stewid.github.io/SimInf/reference/trajectory-SimInf_model-method.md)
for extracting the simulated data, and
[`prevalence`](http://stewid.github.io/SimInf/reference/prevalence-SimInf_model-method.md)
for calculating prevalence from the trajectory.

## Examples

``` r
## For reproducibility, specify the number of threads to use.
set_num_threads(1)

## Create an 'SIR' model with 10 nodes.
model <- SIR(
  u0 = data.frame(
    S = rep(99, 10),
    I = rep(1, 10),
    R = rep(0, 10)
  ),
  tspan = 1:100,
  beta = 0.16,
  gamma = 0.077
)

## Run the model with a fixed seed for reproducibility.
result <- run(model, seed = 22)

## Plot the distribution of susceptible, infected and recovered
## individuals.
plot(result)
```
