# Add spatial coupling information to local data

A utility function to augment local model parameters (`ldata`) with
spatial coupling data from neighboring nodes.

## Usage

``` r
add_spatial_coupling_to_ldata(
  x,
  y,
  cutoff,
  ldata = NULL,
  min_dist = NULL,
  na_fail = TRUE
)
```

## Arguments

- x:

  Numeric vector of projected x coordinates for each node.

- y:

  Numeric vector of projected y coordinates for each node.

- cutoff:

  Numeric scalar. The maximum distance for considering two nodes as
  neighbors. Pairs of nodes farther apart than this value are excluded
  from the neighbor data.

- ldata:

  local data for the nodes. Can either be specified as a `data.frame`
  with one row per node. Or as a matrix where each column `ldata[, j]`
  contains the local data vector for the node `j`. The local data vector
  is passed as an argument to the transition rate functions and the post
  time step function.

- min_dist:

  Numeric scalar. The value to use for the distance between two nodes if
  their coordinates are identical (distance = 0). This prevents division
  by zero errors in downstream calculations (e.g., inverse distance
  weighting). If `NULL` (default) and identical coordinates are found,
  an error is raised.

- na_fail:

  Logical. If `TRUE` (default), missing values (`NA`) in `x` or `y` will
  raise an error. If `FALSE`, distances involving missing coordinates
  are set to zero.

## Value

A numeric matrix with the same number of rows as `ldata`, but with
additional columns containing the neighbor indices and distances as
described in the ‘Output Format’ section.

## Details

The function calculates distances between nodes based on projected
coordinates (`x`, `y`) and appends neighbor data to the `ldata` matrix
for each node.

**Output Format:** The returned matrix has the same number of rows as
the input `ldata`. The columns are organized as follows:

- **Local Parameters**: The first \\n\\ columns correspond to the
  original local parameters passed in `ldata`.

- **Neighbor Pairs**: Following the local parameters, the data is stored
  as pairs of columns: `(neighbor_index, distance)`. The
  `neighbor_index` is a zero-based index of the neighbor node. The
  `distance` is the Euclidean distance to that neighbor.

- **Stop Marker**: Each node's neighbor list is terminated by a pair
  `(-1, 0)` in the position where the next `neighbor_index` would
  appear. This marker appears exactly once per node, immediately after
  the last neighbor.

## See also

[`distance_matrix`](http://stewid.github.io/SimInf/reference/distance_matrix.md)
for computing pairwise distances between nodes, and
[`edge_properties_to_matrix`](http://stewid.github.io/SimInf/reference/edge_properties_to_matrix.md)
for a similar utility that converts edge properties to a matrix format.

## Examples

``` r
## Generate a 5 x 5 grid of nodes separated by 1000m.
nodes <- expand.grid(
  x = seq(from = 0, to = 4000, by = 1000),
  y = seq(from = 0, to = 4000, by = 1000)
)

## Create local data with one parameter per node.
ldata <- matrix(0.1, nrow = 1, ncol = nrow(nodes))

## Add spatial coupling with a 2500m cutoff.
ldata_augmented <- add_spatial_coupling_to_ldata(
  x = nodes$x,
  y = nodes$y,
  cutoff = 2500,
  ldata = ldata
)

## Inspect the result for the first node.
ldata_augmented[, 1]
#>  [1]    0.100    1.000 1000.000    2.000 2000.000    5.000 1000.000    6.000
#>  [9] 1414.214    7.000 2236.068   10.000 2000.000   11.000 2236.068   -1.000
#> [17]    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000
#> [25]    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000
#> [33]    0.000    0.000    0.000    0.000    0.000    0.000    0.000    0.000
#> [41]    0.000    0.000    0.000
```
