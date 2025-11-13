# Display the PMCMC posterior distribution

Display the (approximate) posterior distributions obtained from fitting
a particle Markov chain Monte Carlo algorithm, or the corresponding
trace plots.

## Usage

``` r
# S4 method for class 'SimInf_pmcmc'
plot(x, y, start = 1, end = NULL, thin = 1, ...)
```

## Arguments

- x:

  The `SimInf_pmcmc` object to plot.

- y:

  The trace of all variables and logPost are plotted when `y = "trace"`
  or `y = ~trace`, else the posterior distributions are plotted. Default
  is to plot the posterier distributions.

- start:

  The start iteration to remove some burn-in iterations. Default is
  `start = 1`.

- end:

  the last iteration to include. Default is `NULL` which set `end` to
  the last iteration in the chain.

- thin:

  keep every `thin` iteration after the `start` iteration. Default is
  `thin = 1`, i.e., keep every iteration.

- ...:

  Additional arguments affecting the plot.
