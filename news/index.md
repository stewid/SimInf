# Changelog

## SimInf (development version)

### CHANGES OR IMPROVEMENTS

## SimInf 10.0.0 (2025-11-13)

### BREAKING CHANGES

Backwards incompatible changes that are the reason why the major version
has been incremented.

- Redesigned the S4 class SimInf_pfilter and the interface to using the
  bootstrap filtering algorithm, see the documentation for the `pfilter`
  function. Moreover, a variant of the split-step solver (ssm) was added
  to efficiently simulate multiple particles in one trajectory (mssm).

- The slot `replicates` was added to the `SimInf_model`. The slot
  `replicates` holds the number of replicates of each node in a model.
  This new functionality is used by the bootstrap filtering algorithm.

### CHANGES OR IMPROVEMENTS

- Added functionality to fit models to time series data using the
  Particle Markov Chain Monte Carlo (‘PMCMC’) algorithm of Andrieu and
  others (2010)
  [doi:10.1111/j.1467-9868.2009.00736.x](https://doi.org/10.1111/j.1467-9868.2009.00736.x).

- The model parser (`mparse`) now allows the global data vector `gdata`
  to have parameters without a name.

- On Windows, use gsl from the system, via pkg-config when available.

- Renamed the `abc` function argument `ninit` to `n_init`.

- Renamed the `abc` function argument `npart` to `n_particles`.

- When using `use_enum=TRUE` in `mparse`, the enumeration constants
  `N_COMPARTMENTS_U` and `N_COMPARTMENTS_V` will be automatically added
  to facilitate indexing ‘u’ and ‘v’ in the C code. Additionally,
  enumeration values are added to all enumeration constants.

- The `plot` function of a `SimInf_model` object now accepts a `log`
  argument, see the documentation. Thanks to Alfredo Acosta in PR
  [\#56](https://github.com/stewid/SimInf/issues/56).

- Static code analysis of the C codebase has been performed using the
  `clang-tidy` and `cppcheck` tools in order to improve code.

- The `GNU Indent` program has been used to style the C code for
  consistency and readability.

- Add a function to calculate the Lambert W0 function using the GNU
  Scientific Library (GSL).

## SimInf 9.8.1 (2024-06-21)

CRAN release: 2024-06-21

### CHANGES OR IMPROVEMENTS

- Added package anchors in the documentation for all Rd targets not in
  the package itself and the base packages.

### BUG FIXES

- Fix model parser (mparse) to resolve dependencies between variables in
  the transitions when multiple variables have zero dependencies.

## SimInf 9.7.0 (2024-04-23)

CRAN release: 2024-04-23

### CHANGES OR IMPROVEMENTS

- Two new features have been added to the model parser. First, it is now
  possible to define variables in the transitions, which can then be
  used in calculations of propensities or in calculations of other
  variables. Secondly, enumeration constants can now be generated for
  the indices to each parameter in the `u`, `v`, `ldata`, and `gdata`
  vectors in the generated C code. See the
  [`mparse()`](http://stewid.github.io/SimInf/reference/mparse.md)
  documentation.

- Added a utility function
  [`node_events()`](http://stewid.github.io/SimInf/reference/node_events.md)
  to facilitate cleaning raw individual event data and prepare it for
  usage in SimInf.

- Added a new vignette about how to post-process data from a simulated
  trajectory.

- Added a new section about varying probability of picking individuals
  to the scheduled events vignette.

- Added a missing PROTECT in the C function SimInf_distance_matrix. This
  was uncovered by rchk on CRAN.

## SimInf 9.6.0 (2023-12-20)

CRAN release: 2023-12-20

### CHANGES OR IMPROVEMENTS

- Fixed CRAN check warnings for print format specifiers.

- Added a utility function ‘individual_events’ to facilitate cleaning
  raw individual event data and prepare it for usage in SimInf.

## SimInf 9.5.0 (2023-01-23)

CRAN release: 2023-01-23

### CHANGES OR IMPROVEMENTS

- Fixed the configuration script to use R to find the compiler to use.

## SimInf 9.4.0 (2023-01-06)

CRAN release: 2023-01-06

### CHANGES OR IMPROVEMENTS

- Fix the ‘package_skeleton’ function to generate internal C code with
  valid C entry names if the package name contains ‘.’, for example, for
  a package named ‘pkg.name’.

- Changed the usage of ‘any(is.na(x))’ to ‘anyNA(x)’ in the R code.

- Internal refactoring of the ‘distance_matrix’ function to reduce
  memory usage.

- Ensure to check for a valid model object after updating model data.

- Moved the enumeration of event types to the header file
  ‘inst/SimInf.h’

- Fixed problems identified with static analysis of the C code using the
  cppcheck and scan-build tools.

- Added the getter function ‘u0’ to get the initial compartment state of
  a model.

## SimInf 9.3.1 (2022-10-07)

CRAN release: 2022-10-07

### CHANGES OR IMPROVEMENTS

- Cast the return value from R_GetCCallable(“SimInf”, “SimInf_run”); to
  the correct function pointer type.

- Fix strict-prototype warnings.

## SimInf 9.2.0 (2022-09-03)

CRAN release: 2022-09-03

### CHANGES OR IMPROVEMENTS

- Renamed the set/get funtions `update_u0` and `update_v0` to `u0`and
  `v0`, respectively.

- Some internal refactoring of the C code to set number of threads when
  using OpenMP.

- Changed to coerce to a sparse matrix via virtual classes in order to
  work with the upcoming Matrix package update.

## SimInf 9.1.0 (2022-06-08)

CRAN release: 2022-06-08

### CHANGES OR IMPROVEMENTS

- Added a new built-in Susceptible-Infected-Susceptible model: `SIS`.

- Added the functions `update_u0` and `update_v0` to update the initial
  state in a model.

- Added the `pfilter` function to run a bootstrap particle algorithm on
  a model. See the documentation for an example.

### BUG FIXES

- Fixed the unnamed non-character argument ‘usr’ in par() when ploting
  the density of the ABC posterior distribution.

## SimInf 9.0.0 (2022-04-20)

CRAN release: 2022-04-20

This release of SimInf focuses primarily on improving the functionality
for performing Approximate Bayesian computation.

### BREAKING CHANGES

Backwards incompatible changes that are the reason why the major version
has been incremented.

- Redesigned the S4 class SimInf_abc and the interface to using the
  Approximate Bayesian Computation Sequential Monte Carlo (‘ABC-SMC’)
  algorithm, see the documentation for the ‘abc’ function. Moreover,
  added functionality to adaptively select a sequence of tolerances
  using the algorithm ‘Adaptive Approximate Bayesian Computation
  Tolerance Selection’ of Simola and others (2021), Bayesian Analysis.

### CHANGES OR IMPROVEMENTS

- Added the function ‘n_generations’ to determine the number of
  generations in a ‘SimInf_abc’ object.

### BUG FIXES

- The `plot` function for the SimInf_abc class now passes the additional
  arguments in ‘…’ to the underlying plot functions.

## SimInf 8.4.0 (2021-09-19)

CRAN release: 2021-09-19

### CHANGES OR IMPROVEMENTS

- The ‘events’, ‘gdata’, ‘gdata\<-’, ‘ldata’, ‘punchcard\<-’,
  ‘select_matrix’, ‘select_matrix\<-’ ‘shift_matrix’, and
  ‘shift_matrix\<-’ functions were changed to S4 methods.

### BUG FIXES

- Protect against an integer overflow that could occur in the
  ‘punchcard’ method for a model with many nodes, compartments and
  time-points.

## SimInf 8.3.2 (2021-06-29)

CRAN release: 2021-06-30

### BUG FIXES

- Updated the build configuration script (‘src/Makevars.ucrt’) for
  Windows UCRT to fix the installation failure on CRAN introduced in
  version 8.3.0. Thanks to Tomas Kalibera for providing the patch.

- Improved the test logic for checks using multiple threads to also
  check that the number of available threads is greater than one before
  running a test using two threads. Thanks to Tomas Kalibera for
  reporting a test failure on CRAN r-devel-windows-x86_64-gcc10-UCRT on
  a machine using one thread.

## SimInf 8.3.0 (2021-06-24)

CRAN release: 2021-06-25

### CHANGES OR IMPROVEMENTS

- Changed the R dependency to R(\>= 4.0).

- Refactoring of the build configuration script on Windows.

- Added functionality to fit models to time series data using the
  Approximate Bayesian Computation Sequential Monte Carlo (‘ABC-SMC’)
  algorithm of Toni and others (2009)
  [doi:10.1098/rsif.2008.0172](https://doi.org/10.1098/rsif.2008.0172).

- Added a vignette about scheduled events. This vignette is
  work-in-progress and not yet complete.

## SimInf 8.2.0 (2020-12-06)

CRAN release: 2020-12-06

### IMPROVEMENTS

- Prevalence function: better handling of condition in prevalence
  formula.

- Use ‘Date’ on the x-axis when plotting events if time was specified in
  a Date format.

- Return event ‘time’ as Date and ‘event’ type as character when
  coercing events to a ‘data.frame’ if they were specified as Date and
  character.

## SimInf 8.1.0 (2020-10-18)

CRAN release: 2020-10-18

This release of SimInf focuses primarily on improving the plotting of
trajectory data.

### CHANGES OR IMPROVEMENTS

- The plotting functionality has been improved to allow visualisation of
  the prevalence using the same formula notation as the ‘prevalence’
  function, see the documentation. In addition, the box around the plot
  has been removed and the legend has been moved to the top.

- The prevalence function has been updated to always return a matrix
  when format=‘matrix’. Before it returned a vector instead of a ‘1 x
  length(tspan)’ matrix when the level argument was equal to 1 or 2, or
  when level = 3 and the prevalence included one node only.

## SimInf 8.0.0 (2020-09-13)

CRAN release: 2020-09-13

This release of SimInf includes improvements and changes to facilitate
post-processing of trajectory data and future development of the
package.

### IMPROVEMENTS

- Added the ‘compartments’ and ‘index’ arguments to the ‘boxplot’ and
  ‘pairs’ plotting methods to facilitate analysis of a simulated
  trajectory, see the documentation.

- The build configuration script has been improved to identify if OpenMP
  can be used. It is also possible to skip the check for OpenMP if
  –disable-openmp is specified when installing the package.

### BREAKING CHANGES

Backwards incompatible changes that are the reason why the major version
has been incremented.

- The ‘prevalence’ function was changed to an S4 method. Moreover, the
  ‘type’ character argument was renamed to ‘level’ and changed to an
  integer argument. The reason for renaming the ‘type’ argument was to
  prepare for improvements in the plot function in the future to display
  the prevalence but where there is already an argument called ‘type’.
  Additionally, the ‘as.is’ argument was renamed to ‘format’ for
  clarity. Finally, the unused ‘…’ argument was removed. See the
  documentation for examples.

- The ‘trajectory’ function was changed to an S4 method. Furthermore,
  the node index argument was renamed from ‘node’ to ‘index’ to
  facilitate future development where indices can also represent other
  structures in a model. Additionally, the ‘as.is’ argument was renamed
  to ‘format’ for clarity. Finally, the unused ‘…’ argument was removed.
  See the documentation for examples.

- The ‘Nn’ function to determine the number of nodes in a model has been
  replaced with the S4 method ‘n_nodes’.

## SimInf 7.0.1 (2020-06-18)

CRAN release: 2020-06-18

### CHANGES

- Removed the timestamp from the first line of the C code generated by
  ‘mparse’ in order to keep the hash of the C-code constant, even if the
  code is regenerated at a later time with an identical model. This
  avoids re-compiling the C-code for the model, if the source code has
  already been compiled.

### BUG FIXES

- Fixed a problem that prevented re-compiling the model C-code even
  though it had changed.

- Changed to use the ‘R_FINITE’ instead of ‘isfinite’ in the internal C
  function ‘SimInf_print_status’

## SimInf 7.0.0 (2020-05-23)

CRAN release: 2020-05-23

### CHANGES OR IMPROVEMENTS

- The beta and gamma parameters in the built-in SIR model have been
  moved internally from ‘gdata’ to ‘ldata’ to allow node specific
  parameter values.

- The beta, gamma and epsilon parameters in the built-in SEIR model have
  been moved internally from ‘gdata’ to ‘ldata’ to allow node specific
  parameter values.

- In order to reduce the compilation-time when running multiple
  simulations of a model that contains C code, the MD5 hash of the C
  code is used to determine if a model has already been compiled and the
  DLL loaded, and thus the compilation step can be skipped before
  running a trajectory.

- The C functions ‘SimInf_forward_euler_linear_decay’ and
  ‘SimInf_local_spread’ are now available to be called by C code in
  another package.

- The values in the ‘E’ matrix are now used as weights when sampling
  individuals for the exit, internal transfer and external transfer
  events. The individuals are sampled, one by one, without replacement
  from the compartments specified by ‘E\[, select\]’ in such a way that
  the probability that a particular individual is sampled at a given
  draw is proportional to the weight in ‘E\[, select\]’.

- The values in the ‘E’ matrix are now used as weights for enter events.
  If the column ‘E\[, select\]’ contains several non-zero entries, the
  compartment to add an individual to is sampled in such a way that the
  probability is proportional to the weight in ‘E\[, select\]’.

- The scheduled enter events can now use ‘proportion’ when ‘n = 0’, see
  description of breaking changes below. Additionally, the ‘shift’
  feature can also be used to further control in which compartments
  individuals are added.

### BREAKING CHANGES

Backwards incompatible changes that are the reason why the major version
has been incremented.

- Removed the ‘run_outer’ function.

- Removed the unused ‘threads’ argument from the ‘SimInf_run’ function
  (use ‘set_num_threads’ to specify the number of threads). The
  ‘SimInf_run’ function is the C function that a model calls to simulate
  a trajectory. Because of this change, model packages created with a
  previous version of SimInf must be modified/recreated to work with
  this version of SimInf.

- Events with n = 0 utilize the proportion instead to calculate the
  number of individuals affected by the event. The number was previously
  calculated by multiplying the number of individuals in a node by the
  proportion in the event. This made it tricky to use proportion for a
  scheduled event when the proportion was very small or very large and
  the number of individuals in a node was small, since the result was
  that it always rounded to 0 individuals with a small proportion and
  all individuals with a large proportion. This has been replaced with a
  sampling from a binomial distribution to determine the number of
  individuals affected by the event. Thanks to Thomas Rosendal in PR
  [\#28](https://github.com/stewid/SimInf/issues/28).

## SimInf 6.5.1 (2020-04-01)

CRAN release: 2020-04-01

### BUG FIXES

- Fixed a memory access error in the internal C code that was introduced
  in the recently released 6.5.0 version. Detected by the CRAN gcc-UBSAN
  tests (Undefined Behavior Sanitizer).

## SimInf 6.5.0 (2020-03-29)

CRAN release: 2020-03-29

### IMPROVEMENTS

- It’s now possible to have a ‘tspan’ vector of length one to simulate
  over one time-unit only i.e. \[t, t+1).

- The `trajectory` function has been ported to C and parallelized to
  efficiently transform simulated data from a model to a `data.frame`.

- Static code analysis of the codebase has been performed using the
  `lintr` package in order to improve code style, consistency and
  readability.

- It’s now possible to specify what type of plot should be drawn from
  the simulated data, see the documentation.

- Better documentation of the SISe, SISe3, SISe_sp, and SISe3_sp models.

### CHANGES

- The way to specify the number of threads to use in parallelized
  functions has been changed to fix that specifying the number of
  threads should only affect SimInf and not other packages using OpenMP.
  Use `set_num_threads` to specify the number of threads, see
  documentation. It still works to pass the number of threads to the
  [`run()`](http://stewid.github.io/SimInf/reference/run.md) function,
  however, the `threads` argument will be removed from
  [`run()`](http://stewid.github.io/SimInf/reference/run.md) in a future
  release.

- Updated the vignette to use the ‘set_num_threads’ function.

- To avoid cluttering the error message, the name of the internal
  function that generated the error has been removed from the error
  message (use `traceback` to print the call stack of the last uncaught
  error). Additionally, all error messages now ends with a `.` (full
  stop).

- Removed the row and column names from the internal matrices U and V
  because they were not used anywhere in the code.

### BUG FIXES

- Fix storing solution of the state vectors of the first time point
  until after simulated time has passed the the first time point in
  tspan.

## SimInf 6.4.0 (2019-11-12)

CRAN release: 2019-11-11

### CHANGES

- Updated the vignette, the CITATION file and added the DOI for the JSS
  publication to DESCRIPTION/Description.

- Renamed the NEWS file to NEWS.md and changed to use markdown format
  style.

### BUG FIX

- Removed a timestamp from a test to avoid a possible test failure.

## SimInf 6.3.0 (2019-05-26)

CRAN release: 2019-05-26

### IMPROVEMENTS

- The `prevalence` function now supports an additional condition in the
  formula specification to further control the population to include in
  the calculation, see the documentation.

- It is now possible to write the C code snippet of the `pts_fun`
  argument to `mparse` as a text string with line breaks instead of as a
  vector of character lines.

- The summary function of a model generated with the `mparse`
  functionality now includes the propensity of each transition.

### BREAKING CHANGES

- Removed the `U<-` and `V<-` methods. Use the new `punchcard<-` method
  instead, see the documentation.

### BUG FIXES

- Fix storing solution of the continuous state vector v up to time t,
  but not including t, when t has passed the next time in tspan.

- Fix a link-type optimization type mismatch in the internal C function
  `SimInf_ldata_sp` (detected by the CRAN checks).

## SimInf 6.2.0 (2018-11-20)

CRAN release: 2018-11-20

### IMPROVEMENTS

- The `mparse` model builder functionality has been improved to be able
  to generate models that contain node specific local data and
  continuous state variables. It is also possible to pass the C code of
  the post time step function to `mparse`. See documentation.

- The package skeleton method has been updated to handle the the
  improved functionlity of `mparse`.

- Print summary of the local data when the `ldata` matrix in the model
  has row names.

- It’s now possible to pass `gdata`, `ldata` and `v0` as data.frames to
  the `SimInf_model` function.

- Better error message when an invalid rate or a negative state is
  detected during a simulation.

## SimInf 6.1.0 (2018-08-13)

CRAN release: 2018-08-13

### IMPROVEMENTS

- Print better error message when processing scheduleed events to
  facilitate debugging.

- Some fixes to the vignette.

- Updated CITATION file.

### BUG FIXES

- Fixed broken mparse example.

## SimInf 6.0.0 (2018-04-21)

CRAN release: 2018-04-20

### IMPROVEMENTS

- The `print` and `summary` methods have been improved to include a
  summary of the data in the model and the number of individuals in each
  compartment.

- Added the `trajectory` method to facilitate post-processing of
  simulation data. See documentation.

- The event type, when specifying the scheduled events from a
  data.frame, can now be coded either as a numerical value or as a
  character string, see the documentataion for the `SimInf_events`
  method.

- Now SimInf also includes the All Events Method (AEM) solver (Bauer P.,
  Engblom S. (2015) Sensitivity Estimation and Inverse Problems in
  Spatial Stochastic Models of Chemical Kinetics. In: Abdulle A.,
  Deparis S., Kressner D., Nobile F., Picasso M. (eds) Numerical
  Mathematics and Advanced Applications - ENUMATH 2013. Lecture Notes in
  Computational Science and Engineering, vol 103. Springer, Cham. Doi:
  10.1007/978-3-319-10705-9_51). A core feature of the AEM solver is
  that reaction events are carried out in channels which access private
  streams of random numbers, in contrast to the Gillespie direct method
  where one use only one stream for all events. We have added a new
  `solver` argument to the `run` method.

- The `package_skeleton` method have been updated to handle the new
  `solver` argument.

- The event type, when specifying the scheduled events from a
  data.frame, can now be coded either as a numerical value or as a
  character string, see the documentataion for the `SimInf_events`
  method.

- Added the utility function `select_matrix` to get/set the E matrix of
  the events slot of a SimInf_model object.

- Added the utility function `shift_matrix` to get/set the N matrix of
  the events slot of a SimInf_model object.

- Added the utility function `gdata` to get/set parameters in the global
  data vector that is common to all nodes in a SimInf_model object.

- Added the utility function `ldata` to get the local data vector of a
  node in a SimInf_model object.

### BREAKING CHANGES

- We have changed `mparse` method to make it even easier to specify a
  new model for the Siminf framework, see the documentation.

- Removed the following methods: `U`, `V`, `recovered`, `susceptible`,
  and `infected`. Use the new `trajectory` method instead, see the
  documentation.

- Renamed the type option `bnp` to `nop` in the `prevalence` method to
  calculate the proportion of nodes with at least one case.

- Removed the `seed` argument from the `run` method. Use `set.seed`
  instead.

## SimInf 5.1.0 (2017-10-18)

CRAN release: 2017-10-18

### BUG FIXES

- Added missing calls to `PROTECT` to protect newly created R objects in
  C code from the garbage collector.

- Protect against invalid transition index in direct SSA.

## SimInf 5.0.1

### IMPROVEMENTS

- Improvements in the vignette.

## SimInf 5.0.0 (2017-06-13)

CRAN release: 2017-06-13

### NEW FEATURES

- Added a vignette that provide a technical description of the
  framework, how to use a predefined model in SimInf, demonstrate a case
  study, and finally show how to extend the framework with a user
  defined model.

- Added the `pairs` and `boxplot` methods to visualise the number of
  individuals in each compartment.

- Added `spaghetti` argument to the plot method to draw one line for
  each node.

- Added the method `package_skeleton` to create a package skeleton for a
  model depending on SimInf.

### CHANGES

- Renamed `plot` argument `which` to `N`.

### IMPROVEMENTS

- The core solver now checks the return value from the transition rate
  functions and raise an error if the return value is non-finite or less
  than zero.

- Added example data (`nodes`) with a population of 1600 nodes within a
  50 square kilometer region to facilitate illustration of various
  models.

### BUG FIXES

- Added missing call to `PROTECT` to protect newly created matrices (R
  objects) from the garbage collector.

## SimInf 4.0.0 (2017-03-21)

CRAN release: 2017-03-21

### NEW FEATURES

- Added the `V` and `U` getter/setter methods to the `SimInf_model`.

- The header file `SimInf.h` was added to the folder `inst/include` to
  enable `LinkingTo: SimInf`. A call to `R_RegisterCCallable` was added
  to make the `SimInf_run` C routine available to other packages.

- The slots `U_sparse` and `V_sparse` was added to the `SimInf_model`.
  The slots `U_sparse` and `V_sparse` holds the result if the model was
  initialized to write the solution to a sparse matrix.

- The slot `C_code` was added to the `SimInf_model`. Character vector
  with optional model C code. If non-empty, the C code is written to a
  temporary C-file when the `run` method is called. The temporary C-file
  is compiled and the resulting DLL is dynamically loaded. The DLL is
  unloaded and the temporary files are removed after running the model.

- The `mparse` model parser was added. The parser generates the model
  matrices `G` and `S` and C code from a text string specification of
  the model, for example, c(“S -\> k1*S*I -\> I”, “I -\> k2\*I -\> R”)
  expresses a SIR model. The S4 class `SimInf_mparse` was added to
  support this feature.

### IMPROVEMENTS

- Added the model compartment names to the model `u0` and `U` matrices
  to facilitate inspection of data.

- Removed the max node size constraint when sampling individuals for
  scheduled events from more than two compartments (fixes issue
  [\#3](https://github.com/stewid/SimInf/issues/3)).

- The `tspan` vector argument to the `siminf_model` method can now be
  either an `integer` or a `Date` vector. A `Date` vector is internally
  coerced to a numeric vector as days, where `tspan[1]` becomes the day
  of the year of the first year of the tspan argument. The dates are
  added as names to the numeric vector.

- Improved plot methods by adding colors and increase line widths.

- Improved build configuration to register native routines.

### CHANGES

- The S4 class `siminf_model` was renamed to `SimInf_model`.

- Copy `u0` and `v0` to the first column in `U` and `V`, respectively.
  Prior to this change, the state was copied after the first time step.

## SimInf 3.0.0 (2017-01-29)

CRAN release: 2017-01-29

### NEW FEATURES

- Two new built-in models: `SIR` and `SEIR`.

### CHANGES

- Removed the sub-domain argument. Instead use the local data vector
  `ldata` to pass node specific data.

- Removed the unused parameter `epsilon` from the `SISe_sp` and
  `SISe3_sp` models.

- Refactoring of local spread in the `SISe_sp` and `SISe3_sp` models.

- Added a parameter with a pointer to a random number generator to the
  post time step function.

## SimInf 2.0.0 (2016-05-04)

CRAN release: 2016-05-04

### NEW FEATURES

- Two new built-in models: `SISe_sp` and `SISe3_sp` with spatial
  coupling.

- Added the `run_outer` method to run a model on scaled `gdata`
  parameters.

### IMPROVEMENTS

- Added two synthetic datasets `u0_SISe3` and `events_SISe3` to
  demonstrate how to use SimInf.

- Several improvements to the documentation.

- Added summary method to the siminf_model class and the
  scheduled_events class.

### CHANGES

- The slot N in scheduled_events is changed from dgCMatrix to matrix.

- Refactoring of the `show` method for the siminf_model class and the
  scheduled_events class.

- Renamed function prototype `PostTimeStepFun` to `PTSFun` and
  `PropensityFun` to `TRFun`.

- Renamed `init` argument to `u0` in the `SISe` and `SISe3` methods.

- Changed `node`, `dest`, `select` and `shift` to one-based indices in
  `scheduled_events`.

### BUG FIXES

- Changed to use `gsl_rng_uniform_pos` instead of `gsl_rng_uniform` to
  avoid 0 in the direct SSA. If zero was randomly selected and the first
  compartment empty, the simulator could enter a negative state.

## SimInf 1.0.0 (2016-01-08)

CRAN release: 2016-01-08

- First release.
