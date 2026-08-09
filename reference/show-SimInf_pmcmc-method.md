# Brief summary of a `SimInf_pmcmc` object

Display a summary of a particle Markov chain Monte Carlo (PMCMC)
analysis, including the number of iterations, number of particles,
mixing proportion for the adaptive proposal, acceptance ratio, and
posterior quantiles (2.5 standard deviation for each estimated
parameter.

## Usage

``` r
# S4 method for class 'SimInf_pmcmc'
show(object)
```

## Arguments

- object:

  The
  [`SimInf_pmcmc`](http://stewid.github.io/SimInf/reference/SimInf_pmcmc-class.md)
  object to display.

## Value

The `object`, returned invisibly.

## See also

[`pmcmc`](http://stewid.github.io/SimInf/reference/pmcmc.md) for running
a PMCMC analysis,
[`continue_pmcmc`](http://stewid.github.io/SimInf/reference/continue_pmcmc.md)
for continuing an existing PMCMC chain,
[`SimInf_pmcmc`](http://stewid.github.io/SimInf/reference/SimInf_pmcmc-class.md)
for the class definition.
