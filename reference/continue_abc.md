# Run more generations of ABC SMC

Run more generations of ABC SMC

## Usage

``` r
continue_abc(
  object,
  tolerance = NULL,
  data = NULL,
  verbose = getOption("verbose", FALSE),
  post_gen = NULL
)

# S4 method for class 'SimInf_abc'
continue_abc(
  object,
  tolerance = NULL,
  data = NULL,
  verbose = getOption("verbose", FALSE),
  post_gen = NULL
)
```

## Arguments

- object:

  The `SimInf_abc` object to continue from.

- tolerance:

  A numeric matrix (number of summary statistics \\\times\\ number of
  generations) where each column contains the tolerances for a
  generation and each row contains a sequence of gradually decreasing
  tolerances. Can also be a numeric vector if there is only one summary
  statistic. The tolerance determines the number of generations of
  ABC-SMC to run.

- data:

  Optional data to be passed to the `SimInf_abc@fn` function. Default is
  `NULL`.

- verbose:

  prints diagnostic messages when `TRUE`. The default is to retrieve the
  global option `verbose` and use `FALSE` if it is not set.

- post_gen:

  An optional function that, if non-NULL, is applied after each
  completed generation. The function must accept one argument of type
  `SimInf_abc` with the current state of the fitting process. This
  function can be useful to, for example, save and inspect intermediate
  results.

## Value

A `SimInf_abc` object.
