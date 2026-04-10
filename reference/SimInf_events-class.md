# Class `SimInf_events`

Class to hold data for scheduled events that modify the discrete state
of individuals in a node at a pre-defined time `t`.

## Slots

- `E`:

  The **select matrix** (sparse matrix of class
  [`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)).
  Each row corresponds to a model compartment.

  - **Sampling (Exit, Internal/External Transfer):** Non-zero entries in
    a column indicate which compartments individuals are sampled from.
    The values in `E[, select]` act as **weights** for sampling
    individuals without replacement (probability proportional to
    weight).

  - **Targeting (Enter):** Non-zero entries in a column indicate which
    compartments new individuals are added to. The values in
    `E[, select]` act as **weights** for distributing new individuals
    among the target compartments.

- `N`:

  The **shift matrix** (integer matrix). Determines how individuals are
  moved between compartments during *enter*, *internal transfer*, and
  *external transfer* events.

  - Each row corresponds to a source compartment.

  - Each column corresponds to a specific `shift` value.

  - If `q <- shift`, the entry `N[p, q]` defines the **offset** (number
    of rows to move) for individuals sampled from compartment `p`.

  - The destination compartment is calculated as:
    `destination = p + N[p, q]`.

  - Constraint: `1 <= destination <= number of compartments`.

- `event`:

  Integer vector specifying the event type for each row:

  - `0`: *exit* (remove individuals).

  - `1`: *enter* (add individuals).

  - `2`: *internal transfer* (move within node).

  - `3`: *external transfer* (move between nodes).

  Other values are reserved for future use.

- `time`:

  Integer vector specifying the time step when each event occurs.

- `node`:

  Integer vector specifying the **source node** for the event. For
  *external transfer*, this is the node individuals are moved *from*.
  Range: `1 <= node <= number of nodes`.

- `dest`:

  Integer vector specifying the **destination node** for *external
  transfer* events (individuals moved *to*). For other event types, this
  value is ignored (typically set to 0). Range:
  `1 <= dest <= number of nodes`.

- `n`:

  Integer vector specifying the **number of individuals** affected by
  the event. Must be `n >= 0`.

- `proportion`:

  Numeric vector. If `n[i] == 0`, the number of individuals is sampled
  from a binomial distribution using `proportion[i]` and the current
  population size in the selected compartments. Range:
  `0 <= proportion <= 1`.

- `select`:

  Integer vector specifying which **column** of the matrix `E` to use
  for sampling/targeting for each event. The specific individuals are
  chosen based on the non-zero entries in `E[, select[i]]`.

- `shift`:

  Integer vector specifying which **column** of the matrix `N` to use
  for shifting individuals for each event. Unused for *exit* events.

## See also

[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model.md)
for the main model class that holds the events.
[`SIR`](http://stewid.github.io/SimInf/reference/SIR.md),
[`SEIR`](http://stewid.github.io/SimInf/reference/SEIR.md),
[`SIS`](http://stewid.github.io/SimInf/reference/SIS.md),
[`SISe`](http://stewid.github.io/SimInf/reference/SISe.md) for examples
of how events are passed to model constructors.
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model.md)
(constructor) for details on the `E` and `N` matrices used to define
event behavior. [`run`](http://stewid.github.io/SimInf/reference/run.md)
for executing the simulation with scheduled events. Vignette
`"Scheduled events"` for detailed examples of defining and using events.
