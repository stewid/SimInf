# Detailed summary of a `SimInf_pmcmc` object

Display a detailed summary of a particle Markov chain Monte Carlo
(PMCMC) analysis, including the number of iterations, number of
particles, acceptance ratio, model characteristics (name, number of
nodes, transitions), and posterior quantiles (2.5 75 parameter.

## Usage

``` r
# S4 method for class 'SimInf_pmcmc'
summary(object, ...)
```

## Arguments

- object:

  The
  [`SimInf_pmcmc`](http://stewid.github.io/SimInf/reference/SimInf_pmcmc-class.md)
  object to summarize.

- ...:

  Additional arguments affecting the summary produced. Currently
  ignored.

## Value

`NULL`, returned invisibly.

## Details

Compared to
[`show`](http://stewid.github.io/SimInf/reference/show-SimInf_pmcmc-method.md),
the summary additionally displays model characteristics such as the
model name, number of nodes, and transition details.

## See also

[`pmcmc`](http://stewid.github.io/SimInf/reference/pmcmc.md) for running
a PMCMC analysis,
[`continue_pmcmc`](http://stewid.github.io/SimInf/reference/continue_pmcmc.md)
for continuing an existing PMCMC chain,
[`SimInf_pmcmc`](http://stewid.github.io/SimInf/reference/SimInf_pmcmc-class.md)
for the class definition,
[`show`](http://stewid.github.io/SimInf/reference/show-SimInf_pmcmc-method.md)
for a brief summary.
