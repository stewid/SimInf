# Extract Log-Likelihood

Extract the estimated log-likelihood from a `SimInf_pfilter` object.

## Usage

``` r
# S4 method for class 'SimInf_pfilter'
logLik(object)
```

## Arguments

- object:

  The
  [`SimInf_pfilter`](http://stewid.github.io/SimInf/reference/SimInf_pfilter-class.md)
  object.

## Value

The estimated log-likelihood, returned as a numeric scalar.

## See also

[`pfilter`](http://stewid.github.io/SimInf/reference/pfilter.md) for
running a particle filter and creating `SimInf_pfilter` objects,
[`SimInf_pfilter`](http://stewid.github.io/SimInf/reference/SimInf_pfilter-class.md)
for the class definition.
