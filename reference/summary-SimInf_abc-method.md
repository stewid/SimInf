# Detailed summary of a `SimInf_abc` object

Display a detailed summary of an approximate Bayesian computation (ABC)
analysis, including the number of particles per generation, the number
of generations, and—for each generation—the acceptance rate, effective
sample size (ESS), and summary statistics (minimum, first quartile,
median, mean, third quartile, and maximum) for each parameter in the
posterior distribution.

## Usage

``` r
# S4 method for class 'SimInf_abc'
summary(object, ...)
```

## Arguments

- object:

  The
  [`SimInf_abc`](http://stewid.github.io/SimInf/reference/SimInf_abc-class.md)
  object to summarize.

- ...:

  Additional arguments affecting the summary produced. Currently
  ignored.

## Value

`NULL`, returned invisibly.

## Details

Compared to
[`show`](http://stewid.github.io/SimInf/reference/show-SimInf_abc-method.md),
the summary additionally displays posterior statistics for all
generations, not just the final one.

## See also

[`abc`](http://stewid.github.io/SimInf/reference/abc.md) for running an
ABC analysis,
[`continue_abc`](http://stewid.github.io/SimInf/reference/continue_abc.md)
for continuing an existing ABC chain,
[`SimInf_abc`](http://stewid.github.io/SimInf/reference/SimInf_abc-class.md)
for the class definition,
[`show`](http://stewid.github.io/SimInf/reference/show-SimInf_abc-method.md)
for a brief summary.
