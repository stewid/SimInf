# Detailed summary of a `SimInf_pfilter` object

Display a detailed summary of a bootstrap particle filter analysis,
including the number of particles, log-likelihood estimate, and model
characteristics (model name, number of nodes, scheduled events,
transitions, global data, local data, continuous state variables, and
compartments).

## Usage

``` r
# S4 method for class 'SimInf_pfilter'
summary(object, ...)
```

## Arguments

- object:

  The
  [`SimInf_pfilter`](http://stewid.github.io/SimInf/reference/SimInf_pfilter-class.md)
  object to summarize.

- ...:

  Additional arguments affecting the summary produced. Currently
  ignored.

## Value

`NULL`, returned invisibly.

## Details

Compared to
[`show`](http://stewid.github.io/SimInf/reference/show-SimInf_pfilter-method.md),
the summary additionally displays full model characteristics from the
underlying
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
object.

## See also

[`pfilter`](http://stewid.github.io/SimInf/reference/pfilter.md) for
running a particle filter,
[`logLik`](https://rdrr.io/r/stats/logLik.html) for extracting the
log-likelihood,
[`SimInf_pfilter`](http://stewid.github.io/SimInf/reference/SimInf_pfilter-class.md)
for the class definition,
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model-class.md)
for the model class embedded in the filter,
[`show`](http://stewid.github.io/SimInf/reference/show-SimInf_pfilter-method.md)
for a brief summary.
