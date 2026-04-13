# Model parser to define new models for `SimInf`

Describe your model using a simple string syntax in R. `mparse` parses
this description, generates model-specific C code, and returns a
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object ready for simulation.

## Usage

``` r
mparse(
  transitions = NULL,
  compartments = NULL,
  ldata = NULL,
  gdata = NULL,
  u0 = NULL,
  v0 = NULL,
  tspan = NULL,
  events = NULL,
  E = NULL,
  N = NULL,
  pts_fun = NULL,
  use_enum = FALSE
)
```

## Arguments

- transitions:

  Character vector defining the state transitions. Each element follows
  the format `"Source -> Propensity -> Target"`.

  - **Syntax**: `"X -> rate_expr -> Y"` moves individuals from
    compartment `X` to `Y` with rate `rate_expr`.

  - **Empty Set**: Use `@` for the empty set (e.g., `"I -> mu*I -> @"`
    for death, or `"@ -> lambda -> S"` for birth).

  - **Variables**: Define intermediate variables using `<-`. Example:
    `"N <- S + I + R"`. Variables can be used in propensity expressions.
    Order does not matter.

  - **Types**: By default, variables are `double`. Use `(int)` to force
    integer type (e.g., `"(int)N <- S+I+R"`).

  Example:
  `c("S -> beta*S*I/N -> I", "I -> gamma*I -> R", "N <- S+I+R")`.

- compartments:

  Character vector of compartment names (e.g., `c("S", "I", "R")`).

- ldata:

  Optional local data (node-specific parameters). Can be:

  - A `data.frame` with one row per node (columns = parameters).

  - A numeric matrix where columns are nodes and row names are
    parameters.

  - A named vector (for single-node models).

- gdata:

  Optional global data (common to all nodes). Can be:

  - A named numeric vector (names identify parameters).

  - A one-row `data.frame`.

  Unnamed parameters in a vector must be referenced by index in the
  transitions.

- u0:

  Initial state vector. Can be a `data.frame`, matrix, or named vector.
  See
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  for details.

- v0:

  optional data with the initial continuous state in each node. `v0` can
  be specified as a `data.frame` with one row per node, as a numeric
  matrix where column `v0[, j]` contains the initial state vector for
  the node `j`, or as a named vector when the model only contains one
  node. If `v0` is specified as a `data.frame`, each column is one
  parameter. If `v0` is specified as a matrix, the row names identify
  the parameters. If `v0` is specified as a named vector, the names
  identify the parameters. The ‘v’ vector is passed as an argument to
  the transition rate functions and the post time step function. The
  continuous state can be updated in the post time step function.

- tspan:

  A vector (length \>= 1) of increasing time points where the state of
  each node is to be returned. Can be either an `integer` or a `Date`
  vector.

  - If `integer`: Represents the specific time points (e.g., days,
    hours) at which to record the state.

  - If `Date`: Coerced to a numeric vector representing the **day of the
    year** (1–366) relative to the first date in the vector. The
    original `Date` objects are preserved as names for the numeric
    vector, facilitating time-series plotting.

- events:

  A `data.frame` with the scheduled events. Default is `NULL` i.e. no
  scheduled events in the model.

- E:

  Optional select matrix for events. Can be a `matrix` or `data.frame`.
  Defines which compartments are affected by events and their sampling
  weights. See
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
  for detailed logic on `data.frame` columns (`compartment`, `select`,
  `value`). Default is `NULL` (no events).

- N:

  Optional shift matrix for events. Can be a `matrix` or `data.frame`.
  Defines how individuals are moved between compartments during events.
  See
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
  for detailed logic on `data.frame` columns (`compartment`, `shift`,
  `value`). Default is `NULL` (no events).

- pts_fun:

  optional character vector with C code for the post time step function.
  The C code should contain only the body of the function i.e. the code
  between the opening and closing curly brackets.

- use_enum:

  Logical. If `TRUE`, generates C enumeration constants for parameters
  found in the `u`, `v`, `gdata`, and `ldata` vectors to facilitate
  manual C code modification or use in `pts_fun`. The name of each
  enumeration constant is the upper-case version of the corresponding
  parameter name (e.g., `beta` becomes `BETA`). Additionally, the
  following constants are automatically generated:

  - `N_COMPARTMENTS_U`: The number of discrete compartments.

  - `N_COMPARTMENTS_V`: The number of continuous state variables.

  These constants facilitate indexing of the `u` and `v` vectors in C
  code. Note that `N_COMPARTMENTS_U` and `N_COMPARTMENTS_V` cannot be
  used as compartment or variable names. Default is `FALSE` (parameters
  accessed by integer indices).

## Value

a
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object

## Examples

``` r
## For reproducibility, set the seed.
set.seed(22)

## Create an SIR model with a defined population size variable.
model <- mparse(
  transitions = c(
    "S -> beta * S * I / N -> I",
    "I -> gamma * I -> R",
    "N <- S + I + R"
  ),
  compartments = c("S", "I", "R"),
  gdata = c(beta = 0.16, gamma = 0.077),
  u0 = data.frame(S = 100, I = 1, R = 0),
  tspan = 1:100
)

## Run and plot the result.
result <- run(model)
plot(result)
```
