# SimInf (development version)

## CHANGES

* The name of the internal function that generated the error has been
  removed from the error message (use `traceback` to print the call
  stack of the last uncaught error). Additionally, all error messages
  now ends with a `.` (full stop).

* Removed the row and column names from the internal matrices U and V
  because they were not used anywhere in the code.

## BUG FIXES

* Fix storing solution of the state vectors of the first time point
  until after simulated time has passed the the first time point in
  tspan.

# SimInf 6.3.0

## IMPROVEMENTS

* The `prevalence` function now supports an additional condition in
  the formula specification to further control the population to
  include in the calculation, see the documentation.

* It is now possible to write the C code snippet of the `pts_fun`
  argument to `mparse` as a text string with line breaks instead of as
  a vector of character lines.

* The summary function of a model generated with the `mparse`
  functionality now includes the propensity of each transition.

## BREAKING CHANGES

* Removed the `U<-` and `V<-` methods. Use the new `punchcard<-`
  method instead, see the documentation.

## BUG FIXES

* Fix storing solution of the continuous state vector v up to time t,
  but not including t, when t has passed the next time in
  tspan.

* Fix a link-type optimization type mismatch in the internal C
  function `SimInf_ldata_sp` (detected by the CRAN checks).

# SimInf 6.2.0

## IMPROVEMENTS

* The `mparse` model builder functionality has been improved to be
  able to generate models that contain node specific local data and
  continuous state variables. It is also possible to pass the C code
  of the post time step function to `mparse`. See documentation.

* The package skeleton method has been updated to handle the the
  improved functionlity of `mparse`.

* Print summary of the local data when the `ldata` matrix in the model
  has row names.

* It's now possible to pass `gdata`, `ldata` and `v0` as data.frames
  to the `SimInf_model` function.

* Better error message when an invalid rate or a negative state is
  detected during a simulation.

# SimInf 6.1.0

## IMPROVEMENTS

* Print better error message when processing scheduleed events to
  facilitate debugging.

* Some fixes to the vignette.

* Updated CITATION file.

## BUG FIXES

* Fixed broken mparse example.

# SimInf 6.0.0

## IMPROVEMENTS

* The `print` and `summary` methods have been improved to include a
  summary of the data in the model and the number of individuals in
  each compartment.

* Added the `trajectory` method to facilitate post-processing of
  simulation data. See documentation.

* The event type, when specifying the scheduled events from a
  data.frame, can now be coded either as a numerical value or as a
  character string, see the documentataion for the `SimInf_events`
  method.

* Now SimInf also includes the All Events Method (AEM) solver (Bauer
  P., Engblom S. (2015) Sensitivity Estimation and Inverse Problems in
  Spatial Stochastic Models of Chemical Kinetics. In: Abdulle A.,
  Deparis S., Kressner D., Nobile F., Picasso M. (eds) Numerical
  Mathematics and Advanced Applications - ENUMATH 2013. Lecture Notes
  in Computational Science and Engineering, vol 103. Springer,
  Cham. Doi: 10.1007/978-3-319-10705-9_51).  A core feature of the AEM
  solver is that reaction events are carried out in channels which
  access private streams of random numbers, in contrast to the
  Gillespie direct method where one use only one stream for all
  events. We have added a new `solver` argument to the `run` method.

* The `package_skeleton` method have been updated to handle the new
  `solver` argument.

* The event type, when specifying the scheduled events from a
  data.frame, can now be coded either as a numerical value or as a
  character string, see the documentataion for the `SimInf_events`
  method.

* Added the utility function `select_matrix` to get/set the E matrix
  of the events slot of a SimInf_model object.

* Added the utility function `shift_matrix` to get/set the N matrix
  of the events slot of a SimInf_model object.

* Added the utility function `gdata` to get/set parameters in the
  global data vector that is common to all nodes in a SimInf_model
  object.

* Added the utility function `ldata` to get the local data vector of
  a node in a SimInf_model object.

## BREAKING CHANGES

* We have changed `mparse` method to make it even easier to specify a
  new model for the Siminf framework, see the documentation.

* Removed the following methods: `U`, `V`, `recovered`, `susceptible`,
  and `infected`. Use the new `trajectory` method instead, see the
  documentation.

* Renamed the type option `bnp` to `nop` in the `prevalence` method to
  calculate the proportion of nodes with at least one case.

* Removed the `seed` argument from the `run` method. Use `set.seed`
  instead.

# SimInf 5.1.0

## BUG FIXES

* Added missing calls to `PROTECT` to protect newly created R objects
  in C code from the garbage collector.

* Protect against invalid transition index in direct SSA.

# SimInf 5.0.1

## IMPROVEMENTS

* Improvements in the vignette.

# SimInf 5.0.0

## NEW FEATURES

* Added a vignette that provide a technical description of the
  framework, how to use a predefined model in SimInf, demonstrate a
  case study, and finally show how to extend the framework with a user
  defined model.

* Added the `pairs` and `boxplot` methods to visualise the number of
  individuals in each compartment.

* Added `spaghetti` argument to the plot method to draw one line for
  each node.

* Added the method `package_skeleton` to create a package skeleton for
  a model depending on SimInf.

## CHANGES

* Renamed `plot` argument `which` to `N`.

## IMPROVEMENTS

* The core solver now checks the return value from the transition rate
  functions and raise an error if the return value is non-finite or
  less than zero.

* Added example data (`nodes`) with a population of 1600 nodes within
  a 50 square kilometer region to facilitate illustration of various
  models.

## BUG FIXES

* Added missing call to `PROTECT` to protect newly created matrices (R
  objects) from the garbage collector.

# SimInf 4.0.0

## NEW FEATURES

* Added the `V` and `U` getter/setter methods to the `SimInf_model`.

* The header file `SimInf.h` was added to the folder `inst/include` to
  enable `LinkingTo: SimInf`. A call to `R_RegisterCCallable` was
  added to make the `SimInf_run` C routine available to other
  packages.

* The slots `U_sparse` and `V_sparse` was added to the
  `SimInf_model`. The slots `U_sparse` and `V_sparse` holds the result
  if the model was initialized to write the solution to a sparse
  matrix.

* The slot `C_code` was added to the `SimInf_model`. Character vector
  with optional model C code. If non-empty, the C code is written to a
  temporary C-file when the `run` method is called.  The temporary
  C-file is compiled and the resulting DLL is dynamically loaded. The
  DLL is unloaded and the temporary files are removed after running
  the model.

* The `mparse` model parser was added. The parser generates the model
  matrices `G` and `S` and C code from a text string specification of
  the model, for example, c("S -> k1*S*I -> I", "I -> k2*I -> R")
  expresses a SIR model. The S4 class `SimInf_mparse` was added to
  support this feature.

## IMPROVEMENTS

* Added the model compartment names to the model `u0` and `U` matrices
  to facilitate inspection of data.

* Removed the max node size constraint when sampling individuals for
  scheduled events from more than two compartments (fixes issue #3).

* The `tspan` vector argument to the `siminf_model` method can now be
  either an `integer` or a `Date` vector. A `Date` vector is
  internally coerced to a numeric vector as days, where `tspan[1]`
  becomes the day of the year of the first year of the tspan
  argument. The dates are added as names to the numeric vector.

* Improved plot methods by adding colors and increase line widths.

* Improved build configuration to register native routines.

## CHANGES

* The S4 class `siminf_model` was renamed to `SimInf_model`.

* Copy `u0` and `v0` to the first column in `U` and `V`,
  respectively. Prior to this change, the state was copied after the
  first time step.

# SimInf 3.0.0

## NEW FEATURES

* Two new built-in models: `SIR` and `SEIR`.

## CHANGES

* Removed the sub-domain argument. Instead use the local data vector
  `ldata` to pass node specific data.

* Removed the unused parameter `epsilon` from the `SISe_sp` and
   `SISe3_sp` models.

* Refactoring of local spread in the `SISe_sp` and `SISe3_sp` models.

* Added a parameter with a pointer to a random number generator to the
  post time step function.

# SimInf 2.0.0

## NEW FEATURES

* Two new built-in models: `SISe_sp` and `SISe3_sp` with spatial
  coupling.

* Added the `run_outer` method to run a model on scaled `gdata`
  parameters.

## IMPROVEMENTS

* Added two synthetic datasets `u0_SISe3` and `events_SISe3` to
  demonstrate how to use SimInf.

* Several improvements to the documentation.

* Added summary method to the siminf_model class and the
  scheduled_events class.

## CHANGES

* The slot N in scheduled_events is changed from dgCMatrix to matrix.

* Refactoring of the `show` method for the siminf_model class and the
  scheduled_events class.

* Renamed function prototype `PostTimeStepFun` to `PTSFun` and
  `PropensityFun` to `TRFun`.

* Renamed `init` argument to `u0` in the `SISe` and `SISe3` methods.

* Changed `node`, `dest`, `select` and `shift` to one-based indices in
  `scheduled_events`.

## BUG FIXES

* Changed to use `gsl_rng_uniform_pos` instead of `gsl_rng_uniform` to
  avoid 0 in the direct SSA. If zero was randomly selected and the
  first compartment empty, the simulator could enter a negative state.

# SimInf 1.0.0

* First release.
