# Create a `SimInf_model` object

Construct a low-level `SimInf_model` object. This function is typically
used internally by model constructors (e.g.,
[`SIR()`](http://stewid.github.io/SimInf/reference/SIR.md),
[`mparse()`](http://stewid.github.io/SimInf/reference/mparse.md)) or for
advanced usage where custom model definitions (e.g., user-provided C
code or non-standard matrices) are required.

## Usage

``` r
SimInf_model(
  G,
  S,
  tspan,
  events = NULL,
  ldata = NULL,
  gdata = NULL,
  U = NULL,
  u0 = NULL,
  v0 = NULL,
  V = NULL,
  E = NULL,
  N = NULL,
  C_code = NULL
)
```

## Arguments

- G:

  **Dependency Graph**. Indicates which transition rates need updating
  after a state transition. Can be provided as a sparse matrix (class
  `dgCMatrix`) or a dense matrix. If a dense matrix is provided, it is
  automatically converted to a sparse format internally. See
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  for detailed matrix layout.

- S:

  **State Transition Matrix**. Defines the change in the state vector
  for each transition. Can be provided as a sparse matrix (class
  `dgCMatrix`) or a dense matrix. If a dense matrix is provided, it is
  automatically converted to a sparse format internally. See
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  for detailed matrix layout.

- tspan:

  **Time Span** (numeric or Date vector). Increasing time points for
  output. If `Date`, converted to days with names, where `tspan[1]`
  becomes the day of the year of the first year of `tspan`. The dates
  are added as names to the numeric vector.

- events:

  **Scheduled Events**. A `data.frame` defining the event schedule (see
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)).

- ldata:

  **Local Data**. Parameters specific to each node. Can be:

  - A `data.frame` with one row per node.

  - A matrix where each column `ldata[, j]` is the data vector for node
    `j`.

  Passed to transition rate and post-step functions.

- gdata:

  **Global Data** (numeric vector). Parameters common to all nodes.
  Passed to transition rate and post-step functions.

- U:

  **Result Matrix** (integer matrix). Usually empty at creation. See
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  for detailed matrix layout.

- u0:

  **Initial State**. Initial number of individuals per compartment/node.
  Can be:

  - A matrix (\\N_c \times N_n\\).

  - A `data.frame` with columns corresponding to compartments.

  - Any object coercible to a `data.frame` (e.g., a named numeric vector
    will be coerced to a one-row `data.frame`).

- v0:

  **Initial Continuous State** (numeric matrix). Initial values for
  continuous states per node.

- V:

  **Continuous State Result Matrix** (numeric matrix). Usually empty at
  creation. See
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
  for layout.

- E:

  **Select Matrix** (matrix or `data.frame`). Defines which compartments
  are affected by events and their sampling weights.

  - **Matrix**: Standard sparse matrix.

  - **data.frame**: Must have columns `compartment` and `select`.
    Optional column `value` (default `1`) sets the weight.

  See
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
  for usage details.

- N:

  **Shift Matrix** (matrix or `data.frame`). Defines how individuals are
  moved between compartments during events.

  - **Matrix**: Standard integer matrix.

  - **data.frame**: Must have columns `compartment`, `shift`, and
    `value` (integer offset).

  See
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
  for usage details.

- C_code:

  **C Source Code** (character vector). Optional C code for custom
  transition rates. If provided, it is compiled and loaded when
  [`run()`](http://stewid.github.io/SimInf/reference/run.md) is called.

## Value

A
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object.

## See also

[`SIR`](http://stewid.github.io/SimInf/reference/SIR.md),
[`SEIR`](http://stewid.github.io/SimInf/reference/SEIR.md),
[`SIS`](http://stewid.github.io/SimInf/reference/SIS.md),
[`SISe`](http://stewid.github.io/SimInf/reference/SISe.md) for examples
of compartment model constructors that handle argument validation and
matrix setup.
[`mparse`](http://stewid.github.io/SimInf/reference/mparse.md) for
creating custom models using a simple string syntax.
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for details on the class structure and slots.
[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events.md)
for details on the event schedule format.
