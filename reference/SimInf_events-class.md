# Class `"SimInf_events"`

Class to hold data for scheduled events to modify the discrete state of
individuals in a node at a pre-defined time t.

## Slots

- `E`:

  Each row corresponds to one compartment in the model. The non-zero
  entries in a column indicates the compartments to include in an event.
  For the *exit*, *internal transfer* and *external transfer* events, a
  non-zero entry indicate the compartments to sample individuals from.
  For the *enter* event, all individuals enter first non-zero
  compartment. `E` is sparse matrix of class
  [`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).

- `N`:

  Determines how individuals in *internal transfer* and *external
  transfer* events are shifted to enter another compartment. Each row
  corresponds to one compartment in the model. The values in a column
  are added to the current compartment of sampled individuals to specify
  the destination compartment, for example, a value of `1` in an entry
  means that sampled individuals in this compartment are moved to the
  next compartment. Which column to use for each event is specified by
  the `shift` vector (see below). `N` is an integer matrix.

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

  Determines how individuals in *internal transfer* and *external
  transfer* events are shifted to enter another compartment. The sampled
  individuals are shifted according to column `shift[i]` in matrix `N`
  i.e., `N[, shift[i]]`, where `shift` is an integer vector. See above
  for a description of `N`. Unsued for the other event types.
