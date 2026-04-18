# Length of the MCMC chain

Get the number of iterations (samples) in the Markov Chain Monte Carlo
(MCMC) chain stored in a `SimInf_pmcmc` object.

## Usage

``` r
# S4 method for class 'SimInf_pmcmc'
length(x)
```

## Arguments

- x:

  A `SimInf_pmcmc` object containing the MCMC results.

## Value

An integer scalar representing the number of rows in the `chain` slot
(i.e., the total number of samples in the chain).
