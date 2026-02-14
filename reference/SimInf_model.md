# Create a `SimInf_model`

Create a `SimInf_model`

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

  Dependency graph that indicates the transition rates that need to be
  updated after a given state transition has occured. A non-zero entry
  in element `G[i, i]` indicates that transition rate `i` needs to be
  recalculated if the state transition `j` occurs. Sparse matrix (\\Nt
  \times Nt\\) of object class
  [`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).

- S:

  Each column corresponds to a transition, and execution of state
  transition `j` amounts to adding the `S[, j]` to the state vector of
  the node where the state transition occurred. Sparse matrix (\\Nc
  \times Nt\\) of object class
  [`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).

- tspan:

  A vector (length \>= 1) of increasing time points where the state of
  each node is to be returned. Can be either an `integer` or a `Date`
  vector. A `Date` vector is coerced to a numeric vector as days, where
  `tspan[1]` becomes the day of the year of the first year of `tspan`.
  The dates are added as names to the numeric vector.

- events:

  A `data.frame` with the scheduled events.

- ldata:

  local data for the nodes. Can either be specified as a `data.frame`
  with one row per node. Or as a matrix where each column `ldata[, j]`
  contains the local data vector for the node `j`. The local data vector
  is passed as an argument to the transition rate functions and the post
  time step function.

- gdata:

  A numeric vector with global data that is common to all nodes. The
  global data vector is passed as an argument to the transition rate
  functions and the post time step function.

- U:

  The result matrix with the number of individuals in each disease state
  in every node (\\N_n N_c \times\\ `length(tspan)`). `U[, j]` contains
  the number of individuals in each disease state at `tspan[j]`.
  `U[1:Nc, j]` contains the state of node `1` at `tspan[j]`.
  `U[(Nc + 1):(2 * Nc), j]` contains the state of node `2` at `tspan[j]`
  etc.

- u0:

  The initial state vector. Either a matrix (\\N_c \times N_n\\) or a a
  `data.frame` with the number of individuals in each compartment in
  every node.

- v0:

  The initial continuous state vector in every node. (`dim(ldata)[1]`
  \\\times N_N\\). The continuous state vector is updated by the
  specific model during the simulation in the post time step function.

- V:

  The result matrix for the real-valued continous compartment state
  (\\N_n\\`dim(ldata)[1]` \\\times\\ `length(tspan)`). `V[, j]` contains
  the real-valued state of the system at `tspan[j]`.

- E:

  A matrix to handle scheduled events, see
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md).
  Each row in the matrix corresponds to one compartment in the model.
  The non-zero entries in a column indicates the compartments to include
  in an event. For the *exit*, *internal transfer* and *external
  transfer* events, a non-zero entry indicate the compartments to sample
  individuals from. For the *enter* event, all individuals enter first
  non-zero compartment. The select matrix `E` can either be specified as
  a `matrix`, or as a `data.frame`. When `E` is specified as a
  `data.frame`, it must have one column named `compartment` that defines
  which compartment is referred to, and one column `select` that defines
  the column in `E`. In addition, the `data.frame` can contain an
  optional column named `value` with the value in `E`. When the `value`
  column is missing, `1` is used as the default value.

- N:

  Sparse matrix to handle scheduled events, see
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md).

- C_code:

  Character vector with optional model C code. If non-empty, the C code
  is written to a temporary C-file when the `run` method is called. The
  temporary C-file is compiled and the resulting DLL is dynamically
  loaded. The DLL is unloaded and the temporary files are removed after
  running the model.

## Value

[SimInf_model](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
