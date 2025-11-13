# Lambert W0 function

W0(x) is the principal branch of the solution of the function defined by
\\We^W = x\\ for \\x \>= -1/e\\. The value is calculated using GNU
Scientific Library (GSL).

## Usage

``` r
lambertW0(x)
```

## Arguments

- x:

  numeric vector of values.

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
