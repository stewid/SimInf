# Create a `SISe3` model

Create a `SISe3` model to be used by the simulation framework.

## Usage

``` r
SISe3(
  u0,
  tspan,
  events = NULL,
  phi = NULL,
  upsilon_1 = NULL,
  upsilon_2 = NULL,
  upsilon_3 = NULL,
  gamma_1 = NULL,
  gamma_2 = NULL,
  gamma_3 = NULL,
  alpha = NULL,
  beta_t1 = NULL,
  beta_t2 = NULL,
  beta_t3 = NULL,
  beta_t4 = NULL,
  end_t1 = NULL,
  end_t2 = NULL,
  end_t3 = NULL,
  end_t4 = NULL,
  epsilon = NULL
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

- upsilon_1:

  Indirect transmission rate of the environmental infectious pressure in
  age category 1

- upsilon_2:

  Indirect transmission rate of the environmental infectious pressure in
  age category 2

- upsilon_3:

  Indirect transmission rate of the environmental infectious pressure in
  age category 3

- gamma_1:

  The recovery rate from infected to susceptible for age category 1

- gamma_2:

  The recovery rate from infected to susceptible for age category 2

- gamma_3:

  The recovery rate from infected to susceptible for age category 3

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

- epsilon:

  The background environmental infectious pressure

## Value

`SISe3`

## Details

The `SISe3` model contains two compartments in three age categories;
number of susceptible (S_1, S_2, S_3) and number of infectious (I_1,
I_2, I_3). Additionally, it contains an environmental compartment to
model shedding of a pathogen to the environment. Consequently, the model
has six state transitions,

\$\$S_1 \stackrel{\upsilon_1 \varphi S_1}{\longrightarrow} I_1\$\$

\$\$I_1 \stackrel{\gamma_1 I_1}{\longrightarrow} S_1\$\$

\$\$S_2 \stackrel{\upsilon_2 \varphi S_2}{\longrightarrow} I_2\$\$

\$\$I_2 \stackrel{\gamma_2 I_2}{\longrightarrow} S_2\$\$

\$\$S_3 \stackrel{\upsilon_3 \varphi S_3}{\longrightarrow} I_3\$\$

\$\$I_3 \stackrel{\gamma_3 I_3}{\longrightarrow} S_3\$\$

where the transition rate per unit of time from susceptible to infected
is proportional to the concentration of the environmental contamination
\\\varphi\\ in each node. Moreover, the transition rate from infected to
susceptible is the recovery rate \\\gamma_1, \gamma_2, \gamma_3\\,
measured per individual and per unit of time. Finally, the environmental
infectious pressure in each node is evolved by,

\$\$\frac{d\varphi(t)}{dt} = \frac{\alpha \left(I_1(t) + I_2(t) +
I_3(t)\right)}{N(t)} - \beta(t) \varphi(t) + \epsilon\$\$

where \\\alpha\\ is the average shedding rate of the pathogen to the
environment per infected individual and \\N = S_1 + S_2 + S_3 + I_1 +
I_2 + I_3\\ the size of the node. The seasonal decay and removal of the
pathogen is captured by \\\beta(t)\\. It is also possible to include a
small background infectious pressure \\\epsilon\\ to allow for other
indirect sources of environmental contamination. The environmental
infectious pressure \\\varphi(t)\\ in each node is evolved each time
unit by the Euler forward method. The value of \\\varphi(t)\\ is saved
at the time-points specified in `tspan`.

The argument `u0` must be a `data.frame` with one row for each node with
the following columns:

- S_1:

  The number of sucsceptible in age category 1

- I_1:

  The number of infected in age category 1

- S_2:

  The number of sucsceptible in age category 2

- I_2:

  The number of infected in age category 2

- S_3:

  The number of sucsceptible in age category 3

- I_3:

  The number of infected in age category 3

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
