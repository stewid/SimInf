# Brief summary of a `SimInf_abc` object

Display a summary of an approximate Bayesian computation (ABC) analysis,
including the number of particles per generation, the number of
generations, and—for the final generation—the acceptance rate, effective
sample size (ESS), and summary statistics (minimum, first quartile,
median, mean, third quartile, and maximum) for each parameter in the
posterior distribution.

## Usage

``` r
# S4 method for class 'SimInf_abc'
show(object)
```

## Arguments

- object:

  The
  [`SimInf_abc`](http://stewid.github.io/SimInf/reference/SimInf_abc-class.md)
  object to display.

## Value

The `object`, returned invisibly.

## See also

[`abc`](http://stewid.github.io/SimInf/reference/abc.md) for running an
ABC analysis,
[`continue_abc`](http://stewid.github.io/SimInf/reference/continue_abc.md)
for continuing an existing ABC chain,
[`SimInf_abc`](http://stewid.github.io/SimInf/reference/SimInf_abc-class.md)
for the class definition.
