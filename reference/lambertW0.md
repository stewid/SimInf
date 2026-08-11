# Lambert W0 function

The principal branch of the Lambert W function, which solves \\W e^W =
x\\. Defined for \\x \>= -1/e\\, with output in the range \\\[-1,
\infty)\\. The value is calculated using the GNU Scientific Library
(GSL).

## Usage

``` r
lambertW0(x)
```

## Arguments

- x:

  A numeric vector of values.

## Value

A numeric vector of the same length as `x` containing the principal
branch of the Lambert W function evaluated at each element.

## References

GNU Scientific Library \<https://www.gnu.org/software/gsl/\>

## Examples

``` r
## Should equal 1, as 1 * exp(1) = e.
lambertW0(exp(1))
#> [1] 1

## Should equal 0, as 0 * exp(0) = 0.
lambertW0(0)
#> [1] 0

## Should equal -1.
lambertW0(-exp(-1))
#> [1] -1
```
