# Convert an edge list with properties to a matrix

A utility function to facilitate preparing edge properties for `ldata`
in a model.

## Usage

``` r
edge_properties_to_matrix(edges, n_nodes)
```

## Arguments

- edges:

  a `data.frame` with properties assigned for each edge 'from' –\> 'to',
  for example, weight or count. The `data.frame` must contain the
  columns '`from`' and '`to`' with valid indices to the nodes (1 \<=
  index \<= n_nodes).

- n_nodes:

  the total number of nodes in the model. The resulting matrix will have
  the number of columns equal to `n_nodes`.

## Value

a numeric matrix with the number of rows equal to
`max(table(edges$to)) * (ncol(edges) - 1) + 1` and the number of columns
equal to `n_nodes`.

## Details

The edge properties will be converted to a matrix where each row in
`edges` will become a sequence of (index, value_1, value_2, ...,
value_n) where 'index' is the zero-based index of the `from` node. The
reason for a zero-based index is to facilitate it's usage in C code. The
sequence will be added to the 'to' column in the matrix. There will
always be at least one stop value=-1 in each column. All other values in
the matrix will be set to `NaN`. See ‘Examples’.

## Examples

``` r
## Let us consider the following edge properties.
edges <- data.frame(
    from  = c(  2,    3,     4,  1,   4,    5,   1,   3,   1,   3),
    to    = c(  1,    1,     1,  2,   3,    3,   4,   4,   5,   5),
    rate  = c(0.2, 0.01,  0.79,  1, 0.2, 0.05, 0.2, 0.8, 0.2, 0.8),
    count = c(  5,    5,     5, 50,  10,   10,   5,   5,   5,   5))

## Converting the edge properties into a matrix
edge_properties_to_matrix(edges, 6)
#>        [,1] [,2]  [,3] [,4] [,5] [,6]
#>  [1,]  1.00    0  3.00  0.0  0.0   -1
#>  [2,]  0.20    1  0.20  0.2  0.2  NaN
#>  [3,]  5.00   50 10.00  5.0  5.0  NaN
#>  [4,]  2.00   -1  4.00  2.0  2.0  NaN
#>  [5,]  0.01  NaN  0.05  0.8  0.8  NaN
#>  [6,]  5.00  NaN 10.00  5.0  5.0  NaN
#>  [7,]  3.00  NaN -1.00 -1.0 -1.0  NaN
#>  [8,]  0.79  NaN   NaN  NaN  NaN  NaN
#>  [9,]  5.00  NaN   NaN  NaN  NaN  NaN
#> [10,] -1.00  NaN   NaN  NaN  NaN  NaN

## Gives the following output. The first column contains first the
## properties for the edge from = 2 --> to = 1, where the first
## row is the zero-based index of from, i.e., 1. The second row
## contains the rate=0.2 and the third row count=5. On the fourth
## row starts the next sequence with the values in the second row
## in the edges data.frame. The stop value in the first column is
## on row 10. As can be seen in column 6, there are no edge
## properties for node=6.
##        [,1] [,2]  [,3] [,4] [,5] [,6]
##  [1,]  1.00    0  3.00  0.0  0.0   -1
##  [2,]  0.20    1  0.20  0.2  0.2  NaN
##  [3,]  5.00   50 10.00  5.0  5.0  NaN
##  [4,]  2.00   -1  4.00  2.0  2.0  NaN
##  [5,]  0.01  NaN  0.05  0.8  0.8  NaN
##  [6,]  5.00  NaN 10.00  5.0  5.0  NaN
##  [7,]  3.00  NaN -1.00 -1.0 -1.0  NaN
##  [8,]  0.79  NaN   NaN  NaN  NaN  NaN
##  [9,]  5.00  NaN   NaN  NaN  NaN  NaN
## [10,] -1.00  NaN   NaN  NaN  NaN  NaN
```
