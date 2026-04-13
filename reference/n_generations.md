# Determine the number of generations in an ABC analysis

Extract the number of generations performed in an Approximate Bayesian
Computation (ABC) analysis from a `SimInf_abc` object. Each generation
represents a sequential round of the algorithm, involving the simulation
of particles, acceptance/rejection based on the distance threshold, and
parameter perturbation for the next round.

## Usage

``` r
n_generations(object)

# S4 method for class 'SimInf_abc'
n_generations(object)
```

## Arguments

- object:

  A `SimInf_abc` object containing the results of an ABC analysis.

## Value

An integer scalar representing the total number of generations executed
in the analysis.
