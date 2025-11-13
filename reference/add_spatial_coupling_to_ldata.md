# Add information about spatial coupling between nodes to 'ldata'

A utility function to combine local model parameters (`ldata`) and
spatial coupling to other nodes and add the result to `ldata`.

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

  Projected x coordinate

- y:

  Projected y coordinate

- cutoff:

  The distance cutoff

- ldata:

  local data for the nodes. Can either be specified as a `data.frame`
  with one row per node. Or as a matrix where each column `ldata[, j]`
  contains the local data vector for the node `j`. The local data vector
  is passed as an argument to the transition rate functions and the post
  time step function.

- min_dist:

  The minimum distance to separate two nodes. If the coordinates for two
  nodes are identical, the min_dist must be assigned or an error is
  raised. Default is `NULL`, i.e., to raise an error.

- na_fail:

  A logical indicating whether missing values in `x` or `y` should raise
  an error or assign zero to all distances involving missing values.
  Default is `TRUE`, i.e., to raise an error.

## Value

matrix

## Details

Format for ldata: the first n indicies (1, 2, ..., n) are the local
model parameters, i.e, the indata to the function. They are followed by
the neighbor data, pairs of (index, value) and then a stop pair (-1, 0)
where 'index' is the zero-based index to the neighbor and the value is
determined by the metric argument.

## Examples

``` r
ldata <- add_spatial_coupling_to_ldata(x = nodes$x,
                                       y = nodes$y,
                                       cutoff = 5000)
```
