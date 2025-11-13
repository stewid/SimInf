# Class `"SimInf_model"`

Class to handle data for the `SimInf_model`.

## Slots

- `G`:

  Dependency graph that indicates the transition rates that need to be
  updated after a given state transition has occured. A non-zero entry
  in element `G[i, i]` indicates that transition rate `i` needs to be
  recalculated if the state transition `j` occurs. Sparse matrix (\\Nt
  \times Nt\\) of object class
  [`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).

- `S`:

  Each column corresponds to a state transition, and execution of state
  transition `j` amounts to adding the `S[, j]` column to the state
  vector `u[, i]` of node *i* where the transition occurred. Sparse
  matrix (\\Nc \times Nt\\) of object class
  [`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).

- `U`:

  The result matrix with the number of individuals in each compartment
  in every node. `U[, j]` contains the number of individuals in each
  compartment at `tspan[j]`. `U[1:Nc, j]` contains the number of
  individuals in node 1 at `tspan[j]`. `U[(Nc + 1):(2 * Nc), j]`
  contains the number of individuals in node 2 at `tspan[j]` etc.
  Integer matrix (\\N_n N_c \times\\ `length(tspan)`).

- `U_sparse`:

  If the model was configured to write the solution to a sparse matrix
  ([`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html))
  the `U_sparse` contains the data and `U` is empty. The layout of the
  data in `U_sparse` is identical to `U`. Please note that `U_sparse` is
  numeric and `U` is integer.

- `V`:

  The result matrix for the real-valued continuous state. `V[, j]`
  contains the real-valued state of the system at `tspan[j]`. Numeric
  matrix (\\N_n\\`dim(ldata)[1]` \\\times\\ `length(tspan)`).

- `V_sparse`:

  If the model was configured to write the solution to a sparse matrix
  ([`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html))
  the `V_sparse` contains the data and `V` is empty. The layout of the
  data in `V_sparse` is identical to `V`.

- `ldata`:

  A matrix with local data for the nodes. The column `ldata[, j]`
  contains the local data vector for the node `j`. The local data vector
  is passed as an argument to the transition rate functions and the post
  time step function.

- `gdata`:

  A numeric vector with global data that is common to all nodes. The
  global data vector is passed as an argument to the transition rate
  functions and the post time step function.

- `tspan`:

  A vector of increasing time points where the state of each node is to
  be returned.

- `u0`:

  The initial state vector (\\N_c \times N_n\\) with the number of
  individuals in each compartment in every node.

- `v0`:

  The initial value for the real-valued continuous state. Numeric matrix
  (`dim(ldata)[1]` \\\times N_n\\).

- `events`:

  Scheduled events
  [`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)

- `replicates`:

  Number of replicates of the model.

- `C_code`:

  Character vector with optional model C code. If non-empty, the C code
  is written to a temporary C-file when the `run` method is called. The
  temporary C-file is compiled and the resulting DLL is dynamically
  loaded. The DLL is unloaded and the temporary files are removed after
  running the model.
