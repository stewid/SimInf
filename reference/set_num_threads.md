# Specify the number of threads that SimInf should use

Set the number of threads to be used in SimInf code that is parallelized
with OpenMP (if available). The number of threads is initialized when
SimInf is first loaded in the R session using optional envioronment
variables (see ‘Details’). It is also possible to specify the number of
threads by calling `set_num_threads`. If the environment variables that
affect the number of threads change, then `set_num_threads` must be
called again for it to take effect.

## Usage

``` r
set_num_threads(threads = NULL)
```

## Arguments

- threads:

  integer with maximum number of threads to use in functions that are
  parallelized with OpenMP (if available). Default is NULL, i.e. to use
  all available processors and then check for limits in the environment
  varibles (see ‘Details’).

## Value

The previous value is returned (invisible).

## Details

The `omp_get_num_procs()` function is used to determine the number of
processors that are available to the device at the time the routine is
called. The number of threads is then limited by
`omp_get_thread_limit()` and the current values of the environmental
variables (if set)

- `Sys.getenv("OMP_THREAD_LIMIT")`

- `Sys.getenv("OMP_NUM_THREADS")`

- `Sys.getenv("SIMINF_NUM_THREADS")`

Additionally, the maximum number of threads can be controlled by the
`threads` argument, given that its value is not above any of the limits
described above.
