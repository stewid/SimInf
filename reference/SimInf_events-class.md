# Class `"SimInf_events"`

Class to hold data for scheduled events to modify the discrete state of
individuals in a node at a pre-defined time t.

## Slots

- `E`:

  Each row corresponds to one compartment in the model. The non-zero
  entries in a column indicate the compartments to include in an event.
  For the *exit*, *internal transfer* and *external transfer* events, a
  non-zero entry indicates the compartments to sample individuals from,
  where the values in `E[, select]` are used as weights. Individuals are
  sampled without replacement with probability proportional to the
  weight in `E[, select]`. For the *enter* event, the values in
  `E[, select]` are used as weights when determining which compartment
  to add individuals to. If the column `E[, select]` contains several
  non-zero entries, the compartment is sampled with probability
  proportional to the weight in `E[, select]`. `E` is sparse matrix of
  class
  [`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).

- `N`:

  Determines how individuals in *enter*, *internal transfer* and
  *external transfer* events are shifted to enter another compartment.
  Each row corresponds to one compartment in the model. The values in a
  column define how to move sampled individuals before adding them to
  the destination. Let `q <- shift`, then each non-zero entry in
  `N[, q]` defines the number of rows to move sampled individuals from
  that compartment i.e., sampled individuals from compartment `p` are
  moved to compartment `N[p, q] + p`, where
  `1 <= N[p, q] + p <= N_compartments`. Which column to use for each
  event is specified by the `shift` vector (see below). `N` is an
  integer matrix.

- `event`:

  Type of event: 0) *exit*, 1) *enter*, 2) *internal transfer*, and 3)
  *external transfer*. Other values are reserved for future event types
  and not supported by the current solvers. Integer vector.

- `time`:

  Time of when the event occurs i.e., the event is processed when time
  is reached in the simulation. `time` is an integer vector.

- `node`:

  The node that the event operates on. Also the source node for an
  *external transfer* event. Integer vector. 1 \<= `node[i]` \<= Number
  of nodes.

- `dest`:

  The destination node for an *external transfer* event i.e.,
  individuals are moved from `node` to `dest`, where 1 \<= `dest[i]` \<=
  Number of nodes. Set `event = 0` for the other event types. `dest` is
  an integer vector.

- `n`:

  The number of individuals affected by the event. Integer vector.
  n\[i\] \>= 0.

- `proportion`:

  If `n[i]` equals zero, the number of individuals affected by
  `event[i]` is calculated by sampling the number of individuals from a
  binomial distribution using the `proportion[i]` and the number of
  individuals in the compartments. Numeric vector. 0 \<= proportion\[i\]
  \<= 1.

- `select`:

  To process `event[i]`, the compartments affected by the event are
  specified with `select[i]` together with the matrix `E`, where
  `select[i]` determines which column in `E` to use. The specific
  individuals affected by the event are proportionally sampled from the
  compartments corresponding to the non-zero entries in the specified
  column in `E[, select[i]]`, where `select` is an integer vector.

- `shift`:

  Determines how individuals in *enter*, *internal transfer* and
  *external transfer* events are shifted to enter another compartment.
  The sampled individuals are shifted according to column `shift[i]` in
  matrix `N` i.e., `N[, shift[i]]`, where `shift` is an integer vector.
  Unused for *exit* events.
