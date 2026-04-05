# Getting started with \`mparse\`

## Overview

The `mparse` function is the core engine for defining custom stochastic
disease models in SimInf. Instead of writing complex C code manually,
`mparse` allows you to describe your model’s transitions using a simple,
human-readable string syntax in R. The function then parses this
description, generates optimized C code, compiles it, and returns a
`SimInf_model` object ready for simulation.

This approach offers the best of both worlds: the ease of defining
models in R and the computational speed of compiled C code. It is
particularly powerful for models with complex propensity functions,
multiple compartments, or node-specific parameters.

In this vignette, we will explore:

- The basic syntax for defining transitions.
- How to define variables and population sizes.
- How to incorporate global and local data.
- How to run and visualize the resulting model.

Let us first load the SimInf package.

``` r
library(SimInf)
```

## The Basic Syntax: Transitions

The heart of `mparse` is the `transitions` argument, which is a
character vector describing how individuals move between compartments.
Each transition follows a standard format:

$$\left. \text{Source}\rightarrow\text{Propensity}\rightarrow\text{Destination} \right.$$

- **Source**: The compartment the individual leaves.
- **Propensity**: The rate at which the transition occurs (a
  mathematical expression).
- **Destination**: The compartment the individual enters.

### A Simple SI Model

Let us start with a classic SI model where susceptible individuals ($S$)
become infected ($I$) upon contact. The force of infection is often
modeled as $\beta I/(S + I)$, where $\beta$ is the transmission rate and
$S + I$ is the total population.

We define the transition as follows:

``` r
transitions <- "S -> beta * S * I / (S + I) -> I"
```

Here:

- S is the source compartment.
- beta \* S \* I / (S + I) is the propensity (rate).
- I is the destination compartment.

### Adding Recovery: The SIR Model

To create a standard SIR model, we add a second transition where
infected individuals recover and move to the recovered compartment (R).
The recovery rate is typically $\gamma I$.

``` r
transitions <- c(
  "S -> beta * S * I / (S + I + R) -> I",
  "I -> gamma * I -> R"
)
```

Now we have a complete model definition. Let us create the model object.
We need to specify:

- `transitions`: The vector we just created.
- `compartments`: A vector of all unique compartment names.
- `gdata`: A named vector of global parameters (beta and gamma).
- `u0`: The initial state vector (number of individuals in each
  compartment).
- `tspan`: vector of time points that determines both the duration of
  the simulation and the time points at which the state of the system is
  recorded.

``` r
model <- mparse(
  transitions = transitions,
  compartments = c("S", "I", "R"),
  gdata = c(beta = 0.16, gamma = 0.077),
  u0 = data.frame(S = 99, I = 1, R = 0),
  tspan = 1:100
)
```

The `mparse` function prepares the model definition. Compilation occurs
when [`run()`](http://stewid.github.io/SimInf/reference/run.md) is first
called, and the compiled code is cached for efficiency. Subsequent calls
to [`run()`](http://stewid.github.io/SimInf/reference/run.md) detect the
compiled code and skip the compilation step. Once created, we can run
the simulation and plot the results. For reproducibility, we first call
the [`set.seed()`](https://rdrr.io/r/base/Random.html) function since
there is random sampling involved when picking individuals from the
compartments.

``` r
set.seed(22)
plot(run(model))
```

![Figure 1. Classic SIR epidemic curve generated with
mparse.](mparse_files/figure-html/unnamed-chunk-5-1.png)

Figure 1. Classic SIR epidemic curve generated with mparse.

## Defining Variables and Population Size

In the previous example, we calculated the total population size
directly in the propensity expression as `(S + I + R)`. While this works
perfectly for simple models, repeating the same expression in multiple
transitions can make the code hard to read and maintain. If the
definition of the population changes (e.g., excluding a specific
compartment), you would have to update every transition where it
appears.

`mparse` allows us to define variables using the assignment operator
`<-`. These variables are evaluated within the transition rate functions
and can be reused across multiple transitions. This makes the model
definition cleaner and easier to modify.

### Defining the Total Population

Let us rewrite the SIR model to define the total population `N` as a
variable. We place the variable definition in the `transitions` vector.
The order of definitions does not matter; `mparse` will resolve
dependencies automatically.

``` r
transitions <- c(
  "S -> beta * S * I / N -> I",
  "I -> gamma * I -> R",
  "N <- S + I + R")
```

### Data Types: Integer vs. Double

Although the primary benefit of variables is readability, `mparse` also
allows you to control the data type of the variable. By default,
variables are treated as double (floating-point numbers). However, for
population counts, it is often semantically clearer to define them as
integers. We can enforce this by prefixing the variable name with
`(int)`.

``` r
transitions <- c(
  "S -> beta * S * I / N -> I",
  "I -> gamma * I -> R",
  "(int)N <- S + I + R")
```

Using `(int)N` tells mparse to treat N as an integer. While the compiler
may optimize repeated calculations of `S + I + R` anyway, explicitly
defining N ensures that the logic is centralized and the code remains
easy to read.

### Creating and Running the Model

Let us create the model using the variable definition. Note that we no
longer need to calculate it inline in the propensity expression.

``` r
model <- mparse(
  transitions = transitions,
  compartments = c("S", "I", "R"),
  gdata = c(beta = 0.16, gamma = 0.077),
  u0 = data.frame(S = 99, I = 1, R = 0),
  tspan = 1:100
)
```

Running the model produces the same results as before, but the
transition definitions are now more concise and easier to manage. Note
that we use the same seed value as before.

``` r
set.seed(22)
plot(run(model))
```

![Figure 2. SIR epidemic curve using a defined variable for population
size. The results are identical to Figure
1.](mparse_files/figure-html/unnamed-chunk-9-1.png)

Figure 2. SIR epidemic curve using a defined variable for population
size. The results are identical to Figure 1.

### Handling Edge Cases: Division by Zero

In stochastic simulations, it is possible for a node to become empty
(e.g., all individuals die or move away). If a transition involves
dividing by the total population $N$, and $N$ becomes zero, the
simulation would stop with a “division by zero” error.

Since `mparse` translates the propensity expressions into C code, we can
use the C **ternary operator** (`condition ? true_value : false_value`)
to handle this gracefully.

The syntax `a ? b : c` evaluates to `b` if `a` is true, and `c`
otherwise. For example, to avoid dividing by zero when $N = 0$, we can
write:

``` r
transitions <- c(
  "S -> N > 0 ? beta * S * I / N : 0 -> I",
  "I -> gamma * I -> R",
  "(int)N <- S + I + R"
)
```

Here, the expression `N > 0 ? beta * S * I / N : 0` works as follows:

- If $N > 0$, the force of infection is calculated normally.
- If $N = 0$, the rate is set to 0, preventing the division and allowing
  the simulation to continue (effectively, no new infections can occur
  in an empty population).

This is a robust pattern for any propensity that involves division by a
dynamic quantity, such as the total population size.

### Creating the Model

Let us create the model with this safety check. Note that the logic
remains the same, but the simulation is now robust against empty nodes.

``` r
model <- mparse(
  transitions = transitions,
  compartments = c("S", "I", "R"),
  gdata = c(beta = 0.16, gamma = 0.077),
  u0 = data.frame(S = 99, I = 1, R = 0),
  tspan = 1:100
)
```

``` r
set.seed(22)
plot(run(model))
```

![Figure 3. SIR model using a ternary operator to prevent division by
zero. Although the curve is identical to previous examples (as the
population did not reach zero), this syntax ensures the simulation
continues safely if the node becomes
empty.](mparse_files/figure-html/unnamed-chunk-12-1.png)

Figure 3. SIR model using a ternary operator to prevent division by
zero. Although the curve is identical to previous examples (as the
population did not reach zero), this syntax ensures the simulation
continues safely if the node becomes empty.

## Incorporating Global and Local Data

So far, we have defined models where all parameters (like `beta` and
`gamma`) are the same for the entire population. In many epidemiological
scenarios, however, parameters vary between subpopulations (nodes). For
example, different farms might have different contact rates due to
management practices, while the biological recovery rate remains
constant across all farms.

SimInf distinguishes between two types of data: - **Global Data
(`gdata`)**: Parameters shared by all nodes (e.g., recovery rate
$\gamma$). - **Local Data (`ldata`)**: Parameters specific to each node
(e.g., transmission rate $\beta$).

### A Two-Farm Model

Let us create a model with two independent farms (nodes). We assume: -
Both farms have the same recovery rate (`gamma = 0.077`). - Farm 1 has a
low transmission rate (`beta_farm = 0.1`). - Farm 2 has a high
transmission rate (`beta_farm = 0.4`).

We define the transition using a placeholder name for the local
parameter, `beta_farm`. Note that we also include the safety check for
division by zero using the ternary operator, as discussed previously.

``` r
transitions <- c(
  "S -> N > 0 ? beta_farm * S * I / N : 0 -> I",
  "I -> gamma * I -> R",
  "(int)N <- S + I + R"
)
```

### Defining Global and Local Data

We pass the global parameter `gamma` in gdata. For the local parameter
`beta_farm`, we create a data.frame where each row corresponds to a
node. The column name must match the variable name used in the
transitions.

``` r
gdata <- c(gamma = 0.077)

ldata <- data.frame(
  beta_farm = c(0.1, 0.4)  # Farm 1: 0.1, Farm 2: 0.4
)
```

### Initial Conditions for Multiple Nodes

The initial state `u0` must also be a `data.frame` with one row per
node. The columns correspond to the compartments. The order of rows in
`u0` and `ldata` must match (row 1 is Node 1, row 2 is Node 2, etc.).

It is crucial that the number of rows in `u0` matches the number of rows
in `ldata`. The rows must be aligned by node index: the first row of
both data frames corresponds to Node 1, the second row to Node 2, and so
on. If the row counts differ or the order is mixed, the model will
assign parameters to the wrong nodes, leading to incorrect results.

``` r
u0 <- data.frame(
  S = c(99, 95),  # Farm 1: 99 S, Farm 2: 95 S
  I = c(1, 5),    # Farm 1: 1 I, Farm 2: 5 I
  R = c(0, 0)     # Both start with 0 R
)
```

### Creating and Running the Model

Now we create the model. mparse automatically associates the rows of
ldata and u0 with the nodes (1, 2, …).

``` r
model <- mparse(
  transitions = transitions,
  compartments = c("S", "I", "R"),
  gdata = gdata,
  ldata = ldata,
  u0 = u0,
  tspan = 1:100
)
```

We can now run the simulation. Since we have multiple nodes, the plot()
function can display the trajectories for each node separately if we use
range = FALSE to display the trajectory lines without the shaded range
bands.

``` r
set.seed(22)
result <- run(model)
plot(result, range = FALSE)
```

![Figure 4. Epidemic curves for two farms with different transmission
rates in a single stochastic realization. Farm 2 (higher \`beta_farm\`)
shows a faster outbreak compared to Farm 1, reflecting the expected
impact of the higher transmission rate, though exact timing varies due
to randomness.](mparse_files/figure-html/unnamed-chunk-17-1.png)

Figure 4. Epidemic curves for two farms with different transmission
rates in a single stochastic realization. Farm 2 (higher `beta_farm`)
shows a faster outbreak compared to Farm 1, reflecting the expected
impact of the higher transmission rate, though exact timing varies due
to randomness.

In this realization, Farm 2 (the one with the higher `beta_farm`) shows
a faster outbreak compared to Farm 1. On average, the higher
transmission rate leads to a steeper rise in infections, though
individual stochastic runs may vary. This demonstrates how `ldata`
allows us to easily simulate heterogeneity across a network of
populations without rewriting the transition logic.

### Inspecting Results by Node

To inspect the results for a specific node, we can use the `index`
argument in [`plot()`](https://rdrr.io/r/graphics/plot.default.html) or
[`trajectory()`](http://stewid.github.io/SimInf/reference/trajectory.md).

``` r
plot(result, index = 2)
```

![Figure 5. Trajectory for Farm 2 only, showing the rapid spread due to
the high transmission rate in this specific
realization.](mparse_files/figure-html/unnamed-chunk-18-1.png)

Figure 5. Trajectory for Farm 2 only, showing the rapid spread due to
the high transmission rate in this specific realization.

This flexibility makes mparse ideal for spatial models where each node
has unique characteristics, such as different herd sizes, management
practices, or risk factors.
