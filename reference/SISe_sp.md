# Create a `SISe_sp` model

Create a `SISe_sp` model to be used by the simulation framework.

## Usage

``` r
SISe_sp(
  u0,
  tspan,
  events = NULL,
  phi = NULL,
  upsilon = NULL,
  gamma = NULL,
  alpha = NULL,
  beta_t1 = NULL,
  beta_t2 = NULL,
  beta_t3 = NULL,
  beta_t4 = NULL,
  end_t1 = NULL,
  end_t2 = NULL,
  end_t3 = NULL,
  end_t4 = NULL,
  coupling = NULL,
  distance = NULL
)
```

## Arguments

- u0:

  A `data.frame` with the initial state in each node, i.e., the number
  of individuals in each compartment in each node when the simulation
  starts (see ‘Details’). The parameter `u0` can also be an object that
  can be coerced to a `data.frame`, e.g., a named numeric vector will be
  coerced to a one row `data.frame`.

- tspan:

  A vector (length \>= 1) of increasing time points where the state of
  each node is to be returned. Can be either an `integer` or a `Date`
  vector. A `Date` vector is coerced to a numeric vector as days, where
  `tspan[1]` becomes the day of the year of the first year of `tspan`.
  The dates are added as names to the numeric vector.

- events:

  a `data.frame` with the scheduled events, see
  [`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model.md).

- phi:

  A numeric vector with the initial environmental infectious pressure in
  each node. Will be repeated to the length of nrow(u0). Default is NULL
  which gives 0 in each node.

- upsilon:

  Indirect transmission rate of the environmental infectious pressure

- gamma:

  The recovery rate from infected to susceptible

- alpha:

  Shed rate from infected individuals

- beta_t1:

  The decay of the environmental infectious pressure in interval 1.

- beta_t2:

  The decay of the environmental infectious pressure in interval 2.

- beta_t3:

  The decay of the environmental infectious pressure in interval 3.

- beta_t4:

  The decay of the environmental infectious pressure in interval 4.

- end_t1:

  vector with the non-inclusive day of the year that ends interval 1 in
  each node. Will be repeated to the length of nrow(u0).

- end_t2:

  vector with the non-inclusive day of the year that ends interval 2 in
  each node. Will be repeated to the length of nrow(u0).

- end_t3:

  vector with the non-inclusive day of the year that ends interval 3 in
  each node. Will be repeated to the length of nrow(u0).

- end_t4:

  vector with the non-inclusive day of the year that ends interval 4 in
  each node. Will be repeated to the length of nrow(u0).

- coupling:

  The coupling between neighboring nodes

- distance:

  The distance matrix between neighboring nodes

## Value

`SISe_sp`

## Details

The `SISe_sp` model contains two compartments; number of susceptible (S)
and number of infectious (I). Additionally, it contains an environmental
compartment to model shedding of a pathogen to the environment.
Moreover, it also includes a spatial coupling of the environmental
contamination among proximal nodes to capture between-node spread
unrelated to moving infected individuals. Consequently, the model has
two state transitions,

\$\$S \stackrel{\upsilon \varphi S}{\longrightarrow} I\$\$

\$\$I \stackrel{\gamma I}{\longrightarrow} S\$\$

where the transition rate per unit of time from susceptible to infected
is proportional to the concentration of the environmental contamination
\\\varphi\\ in each node. Moreover, the transition rate from infected to
susceptible is the recovery rate \\\gamma\\, measured per individual and
per unit of time. Finally, the environmental infectious pressure in each
node is evolved by,

\$\$\frac{d \varphi_i(t)}{dt} = \frac{\alpha I\_{i}(t)}{N_i(t)} +
\sum_k{\frac{\varphi_k(t) N_k(t) - \varphi_i(t) N_i(t)}{N_i(t)} \cdot
\frac{D}{d\_{ik}}} - \beta(t) \varphi_i(t)\$\$

where \\\alpha\\ is the average shedding rate of the pathogen to the
environment per infected individual and \\N = S + I\\ the size of the
node. Next comes the spatial coupling among proximal nodes, where \\D\\
is the rate of the local spread and \\d\_{ik}\\ the distance between
holdings \\i\\ and \\k\\. The seasonal decay and removal of the pathogen
is captured by \\\beta(t)\\. The environmental infectious pressure
\\\varphi(t)\\ in each node is evolved each time unit by the Euler
forward method. The value of \\\varphi(t)\\ is saved at the time-points
specified in `tspan`.

The argument `u0` must be a `data.frame` with one row for each node with
the following columns:

- S:

  The number of sucsceptible

- I:

  The number of infected

## Beta

The time dependent beta is divided into four intervals of the year

    where 0 <= day < 365

    Case 1: END_1 < END_2 < END_3 < END_4
    INTERVAL_1 INTERVAL_2     INTERVAL_3     INTERVAL_4     INTERVAL_1
    [0, END_1) [END_1, END_2) [END_2, END_3) [END_3, END_4) [END_4, 365)

    Case 2: END_3 < END_4 < END_1 < END_2
    INTERVAL_3 INTERVAL_4     INTERVAL_1     INTERVAL_2     INTERVAL_3
    [0, END_3) [END_3, END_4) [END_4, END_1) [END_1, END_2) [END_2, 365)

    Case 3: END_4 < END_1 < END_2 < END_3
    INTERVAL_4 INTERVAL_1     INTERVAL_2     INTERVAL_3     INTERVAL_4
    [0, END_4) [END_4, END_1) [END_1, END_2) [END_2, END_3) [END_3, 365)
