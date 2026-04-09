# Specify the number of threads that SimInf should use

Set the number of threads to be used in SimInf code that is parallelized
with OpenMP (if available).

## Usage

``` r
set_num_threads(threads = NULL)
```

## Arguments

- threads:

  Integer specifying the maximum number of threads to use in
  OpenMP-parallelized functions. If `NULL` (default), SimInf attempts to
  use all available processors, subject to the limits imposed by the
  environment variables listed above.

## Value

The previous value of the thread count is returned invisibly.

## Details

The number of threads is initialized when SimInf is first loaded in the
R session, based on optional environment variables (see ‘Details’). It
can also be explicitly set by calling `set_num_threads`. If the
environment variables affecting the thread count change,
`set_num_threads` must be called again for the new values to take
effect.

The function determines the number of available processors using
`omp_get_num_procs()` and limits the thread count based on
`omp_get_thread_limit()` and the following environment variables
(checked in order of precedence):

- `SIMINF_NUM_THREADS`: Specific limit for SimInf.

- `OMP_NUM_THREADS`: General OpenMP limit.

- `OMP_THREAD_LIMIT`: Maximum thread limit.

The `threads` argument allows you to override these limits, provided the
requested value does not exceed the maximum allowed by the environment
or hardware constraints.
