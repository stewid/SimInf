# Create a distance matrix between nodes for spatial models

Calculate the Euclidean distances between all pairs of nodes based on
their projected coordinates (`x`, `y`). Distances greater than the
specified `cutoff` are excluded from the result (stored as zeros in the
sparse matrix).

## Usage

``` r
distance_matrix(x, y, cutoff, min_dist = NULL, na_fail = TRUE)
```

## Arguments

- x:

  Numeric vector of projected x coordinates for each node.

- y:

  Numeric vector of projected y coordinates for each node.

- cutoff:

  Numeric scalar. The maximum distance to include in the matrix. Pairs
  of nodes farther apart than this value are excluded from the sparse
  structure (stored as zeros).

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

A symmetric sparse matrix of class
[`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).
Non-zero entries represent distances \\\le\\ `cutoff`; entries outside
the cutoff are implicitly zero.

## Details

The result is a symmetric sparse matrix (`dgCMatrix`) where the element
`d[i, j]` contains the distance between node `i` and node `j` if it is
less than or equal to `cutoff`, and `0` otherwise.

## Examples

``` r
## Generate a 10 x 10 grid of nodes separated by 100m.
nodes <- expand.grid(
  x = seq(from = 0, to = 900, by = 100),
  y = seq(from = 0, to = 900, by = 100)
)
plot(nodes, main = "Node Grid")


## Calculate distances with a 300m cutoff.
## Only neighbors within 300m will have non-zero entries.
d <- distance_matrix(
  x = nodes$x,
  y = nodes$y,
  cutoff = 300
)

## Inspect the sparse matrix structure.
## Note: The matrix is symmetric and the diagonal is zero.
d[1:10, 1:10]
#> 10 x 10 sparse Matrix of class "dgCMatrix"
#>                                              
#>  [1,]   . 100 200 300   .   .   .   .   .   .
#>  [2,] 100   . 100 200 300   .   .   .   .   .
#>  [3,] 200 100   . 100 200 300   .   .   .   .
#>  [4,] 300 200 100   . 100 200 300   .   .   .
#>  [5,]   . 300 200 100   . 100 200 300   .   .
#>  [6,]   .   . 300 200 100   . 100 200 300   .
#>  [7,]   .   .   . 300 200 100   . 100 200 300
#>  [8,]   .   .   .   . 300 200 100   . 100 200
#>  [9,]   .   .   .   .   . 300 200 100   . 100
#> [10,]   .   .   .   .   .   . 300 200 100   .

## Count the number of neighbors for the first node.
sum(d[1, ] > 0)
#> [1] 10
```
