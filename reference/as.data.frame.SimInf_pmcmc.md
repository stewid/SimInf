# Coerce a `SimInf_pmcmc` object to a `data.frame`

Extract the posterior samples from the MCMC chain stored in a
`SimInf_pmcmc` object and convert them into a `data.frame`.

## Usage

``` r
# S3 method for class 'SimInf_pmcmc'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `SimInf_pmcmc` object.

- ...:

  Additional arguments (currently ignored).

## Value

A `data.frame` where rows represent MCMC iterations and columns
represent the posterior samples of the parameters.

## Details

The resulting `data.frame` contains one row per MCMC iteration and one
column per parameter. These samples represent the joint posterior
distribution of the parameters. This format is convenient for
post-processing and visualization.
