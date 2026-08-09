# Brief summary of a `SimInf_pfilter` object

Display a brief summary of a bootstrap particle filter analysis,
including the number of particles used and the estimated log-likelihood.

## Usage

``` r
# S4 method for class 'SimInf_pfilter'
show(object)
```

## Arguments

- object:

  The
  [`SimInf_pfilter`](http://stewid.github.io/SimInf/reference/SimInf_pfilter-class.md)
  object to display.

## Value

The `object`, returned invisibly.

## See also

[`pfilter`](http://stewid.github.io/SimInf/reference/pfilter.md) for
running a particle filter,
[`logLik`](https://rdrr.io/r/stats/logLik.html) for extracting the
log-likelihood,
[`SimInf_pfilter`](http://stewid.github.io/SimInf/reference/SimInf_pfilter-class.md)
for the class definition.
